#!/usr/bin/env python
#File dedicated to impement the TJPCov pipeline stage into ceci
from .file_types import SACCFile, FiducialCosmology
from ceci import PipelineStage
from ceci.config import StageParameter
import numpy as np

# Some older CROW/TJPCov code paths still reference np.bool (removed in
# recent numpy). Alias once at import time rather than per-call.
if not hasattr(np, 'bool'):
    np.bool = bool


class TJPCovPipeline(PipelineStage):
    """
    TJPCov pipeline stage for covariance computation.

    This stage:
    - Reads an input SACC file (data vector)
    - Reads the fiducial cosmology (shared with TXPipe)
    - Computes covariance terms using TJPCov
    - Optionally merges in non-cluster covariance blocks from the input file
    - Optionally replaces cluster-count covariance using CROW
    - Writes a single output SACC file with the final covariance

    Key configuration groups:
    - Covariance selection (cov_type)
    - Mass–observable relation (mor_parameters)
    - Pipeline behavior (replace_tjpcov_cov)
    """
    name = "TJPCovPipeline"

    inputs = [
        ("clusters_sacc_file", SACCFile),  # For firecrown Likelihood
        ("fiducial_cosmology", FiducialCosmology),
    ]

    outputs = [
        ("clusters_sacc_file_cov", SACCFile),
    ]
    config_options = {
        "replace_tjpcov_cov": StageParameter(bool, True, msg="Replace TJPCov cluster-count covariance with CROW covariance."),
        # Selection function
        "sel_func": StageParameter(bool, True, msg="Include purity and completeness selection functions."),
        "wazp_catalog": StageParameter(bool, False, msg="Use WaZP completeness parameters instead of Aguena+16 defaults (extra, non-default)."),
        "diagonal_shear_covariance": StageParameter(
            bool, True,
            msg=(
                "If True, keep only the diagonal (per radius-bin variance) of "
                "non-cluster-count covariance blocks (e.g. cluster_delta_sigma) "
                "when merging them from the input SACC file."
            ),
        ),
        # TJPCov options
        "use_mpi": StageParameter(bool, False, msg="Use MPI parallelization in TJPCov."),
        "do_xi": StageParameter(bool, False, msg="Compute xi covariance terms."),
        "cov_type": StageParameter(list, ["ClusterCountsGaussian", "ClusterCountsSSC"], msg="TJPCov covariance terms to compute."),
        # Photo-z
        "photo-z": StageParameter(dict, {}, msg="Photo-z uncertainty parameters."),
        # Mass-observable relation
        "mor_parameters": StageParameter(dict, {}, msg="Mass-observable relation parameters."),
    }

    def run(self):
        """
        Main execution:

        - Load input SACC file and fiducial cosmology
        - Compute covariance via TJPCov
        - Merge with existing non-cluster covariance blocks (if needed)
        - Optionally replace cluster-count covariance
        - Save the final result exactly once
        """
        import sacc
        import time
        import os

        # tjpcov packages
        from tjpcov.covariance_calculator import CovarianceCalculator

        st = time.time()
        tjpcov_out_sacc = self.get_output('clusters_sacc_file_cov', final_name=True)
        config_dict = self.config.to_dict()
        outdir = os.path.dirname(tjpcov_out_sacc)
        filename = os.path.basename(tjpcov_out_sacc)
        config_dict['outdir'] = outdir
        sacc_file = self.get_input("clusters_sacc_file")
        # Load the SACC file
        sacc_obj = sacc.Sacc.load_fits(sacc_file)

        with self.open_input("fiducial_cosmology", wrapper=True) as f:
            ccl_cosmo = f.to_ccl()
 
        config_dict['cosmo'] = ccl_cosmo
        ## This duplication below is because of code that should be fixed in TJPCov
        config_dict['parameters'] = {
            "Omega_c": ccl_cosmo["Omega_c"],
            "Omega_b": ccl_cosmo["Omega_b"],
            "h": ccl_cosmo["h"],
            "n_s": ccl_cosmo["n_s"],
            "sigma8": ccl_cosmo["sigma8"],
            "w0": ccl_cosmo["w0"],
            "wa": ccl_cosmo["wa"],
            "Omega_k": ccl_cosmo["Omega_k"],
        }

        # Check if the data contains only cluster counts
        data_types = set(dp.data_type for dp in sacc_obj.data)
        only_counts = data_types == {sacc.standard_types.cluster_counts}

        # Check if covariance is present
        has_covariance = sacc_obj.covariance is not None

        # TJPCov expects settings both nested under "tjpcov" (most of
        # CovarianceBuilder) and at top level (e.g. get_cosmology() reads
        # config["cosmo"] directly).
        config_dict['sacc_file'] = sacc_file
        combined_config = {'tjpcov': config_dict}
        combined_config.update(config_dict)

        cc = CovarianceCalculator(combined_config)
        cov_terms = cc.get_covariance_terms()   # {'gauss': array, 'SSC': array}
        full_cov = cc.get_covariance()
        sacc_with_cov = cc.create_sacc_cov(output=filename, save_terms=True)
        print('Time: ', (time.time() - st), ' sec')

        # From here on, full_cov and sacc_with_cov are the single source of
        # truth. Both post-processing steps mutate and hand them forward;
        # nothing is reloaded from disk, and we save exactly once at the end.
        if not only_counts and has_covariance:
            full_cov = self.merge_data_covariance(
                sacc_obj, full_cov,
                diagonal_only=config_dict.get("diagonal_shear_covariance", True),
            )

        if config_dict["replace_tjpcov_cov"]:
            print("Replacing Counts TJPCov cov for Crow cov. SSC + Counts")
            full_cov = self.replace_crow_counts(
                config_dict, sacc_with_cov, cov_terms, full_cov, ccl_cosmo
            )

        sacc_with_cov.covariance = sacc.covariance.FullCovariance(full_cov)
        sacc_with_cov.save_fits(tjpcov_out_sacc, overwrite=True)

    def extract_and_save_cluster_counts(self, input_sacc_file: str, output_sacc_file: str):
        """
        Reads a SACC file, extracts only the cluster counts data (without covariance),
        and saves it to a new SACC file.

        Args:
            input_sacc_file (str): Path to the input SACC file containing full data.
            output_sacc_file (str): Path where the new SACC file with only cluster counts will be saved.
        """
        import sacc
        # Load the input SACC file
        sacc_obj = sacc.Sacc.load_fits(input_sacc_file)

        # Create a new SACC object for the extracted data
        new_sacc = sacc.Sacc()

        # Copy relevant tracers to the new SACC object
        for tracer_name, tracer in sacc_obj.tracers.items():
            new_sacc.add_tracer_object(tracer)

        # Extract only cluster counts data points
        cluster_count_type = sacc.standard_types.cluster_counts
        cluster_counts_points = sacc_obj.get_data_points(cluster_count_type)

        for point in cluster_counts_points:
            new_sacc.add_data_point(cluster_count_type, point.tracers, point.value)

        # Save the new SACC object to the specified file (without covariance)
        new_sacc.to_canonical_order()
        new_sacc.save_fits(output_sacc_file, overwrite=True)
        return output_sacc_file

    def merge_data_covariance(self, sacc_obj, full_cov, diagonal_only=True):
        """
        Copy non-cluster-count covariance blocks from the original input
        SACC file's covariance into the newly computed full covariance.

        TJPCov only computes cluster-related covariance terms in this
        pipeline (cov_type is restricted to Cluster* classes), so any other
        data type (e.g. cluster_delta_sigma) that was already present with a
        covariance in the input file needs to be preserved here rather than
        left at TJPCov's placeholder value.

        Args:
            sacc_obj (sacc.Sacc): The original input SACC file, still
                carrying its original covariance matrix.
            full_cov (np.ndarray): The newly computed full covariance array
                to merge non-cluster blocks into (modified in place and
                returned).
            diagonal_only (bool): If True, only copy the diagonal (per-point
                variance) for each non-cluster-count data type, dropping any
                off-diagonal correlation between data points of that type.

        Returns:
            np.ndarray: full_cov with non-cluster-count blocks copied over
            from sacc_obj, either diagonal-only or full dense blocks
            depending on diagonal_only.
        """
        import sacc

        data_types_sacc = [
            d_type for d_type in sacc_obj.get_data_types()
            if d_type != sacc.standard_types.cluster_counts
        ]
        for d_type in data_types_sacc:
            ix1 = sacc_obj.indices(data_type=d_type)
            if diagonal_only:
                # Deliberately keep only per-point variance: off-diagonal
                # radius-radius terms are dropped, not accidentally lost.
                full_cov[ix1, ix1] = sacc_obj.covariance.covmat[ix1, ix1]
            else:
                # np.ix_ is required for a full block copy: covmat[ix1, ix1]
                # would only pick out diagonal elements pairwise, not the
                # full (ix1 x ix1) block.
                full_cov[np.ix_(ix1, ix1)] = sacc_obj.covariance.covmat[np.ix_(ix1, ix1)]
        return full_cov

    def replace_crow_counts(self, config_dict, sacc_full, cov_terms, full_cov, cosmo):
        """
        Replace TJPCov cluster-count covariance using CROW predictions.

        This is a temporary workaround.

        Steps:
        - Construct mass–observable relation (mor_parameters)
        - Compute theoretical counts
        - Replace covariance elements using SSC scaling

        WARNING:
        - Hardcoded modeling choices (mass function, grids)
        - Should eventually be implemented inside TJPCov

        Args:
            config_dict (dict): Stage configuration.
            sacc_full (sacc.Sacc): SACC object with tracers/data points for
                the full covariance (used to look up tracer metadata).
            cov_terms (dict): {'gauss': array, 'SSC': array} raw covariance
                term arrays from TJPCov.
            full_cov (np.ndarray): Full covariance array to modify in place.
            cosmo (pyccl.Cosmology): fiducial cosmology, the same object
                used for TJPCov's own covariance computation in run().

        Returns:
            np.ndarray: full_cov with cluster-count blocks replaced.
        """
        import sacc
        import pyccl as ccl
        from crow import ClusterAbundance
        from crow.recipes.binned_grid import GridBinnedClusterRecipe
        from crow import completeness_models, mass_proxy, purity_models, kernel

        # This function should not exist as it should be implemented in TJPCov.
        # This is temporary and so most of the options and configurations are fixed.
        is_wazp = config_dict.get("wazp_catalog", False)
        sel_func = config_dict.get("sel_func", True)
        mor_params = config_dict["mor_parameters"]
        hmf = ccl.halos.MassFuncDespali16(mass_def="200c")
        mass_richness_unb = mass_proxy.MurataUnbinned(pivot_log_mass=mor_params["m_pivot"], pivot_redshift=mor_params["z_pivot"])
        mass_richness_unb.parameters["mu0"]    = mor_params["mu_p0"]
        mass_richness_unb.parameters["mu1"]    = mor_params["mu_p1"]
        mass_richness_unb.parameters["mu2"]    = mor_params["mu_p2"]
        mass_richness_unb.parameters["sigma0"] = mor_params["sigma_p0"]
        mass_richness_unb.parameters["sigma1"] = mor_params["sigma_p1"]
        mass_richness_unb.parameters["sigma2"] = mor_params["sigma_p2"]

        mass_grid_size     = 80
        redshift_grid_size = 40
        proxy_grid_size    = 40
        mass_interval      = (np.log10(float(mor_params["min_halo_mass"])), np.log10(float(mor_params["max_halo_mass"])))
        cl_abundance          = ClusterAbundance(cosmo, hmf)
        purity_aguena         = None#purity_models.PurityAguena16LnProxy()
        completeness_aguena   = completeness_models.CompletenessAguena16()
        if is_wazp:
            purity_aguena = None
            completeness_aguena['a_n'] = 1.570597
            completeness_aguena['b_n'] = -0.028690
            completeness_aguena['a_logm_piv'] = 14.264386
            completeness_aguena['b_logm_piv'] = 0.029814
        redshift_distribution = kernel.SpectroscopicRedshift()
        if not sel_func:
            completeness_aguena = None
            purity_aguena = None
        recipe_grid_abundance = GridBinnedClusterRecipe(
            mass_interval=mass_interval,
            cluster_theory=cl_abundance,
            redshift_distribution=redshift_distribution,
            mass_distribution=mass_richness_unb,
            proxy_grid_size=proxy_grid_size,
            redshift_grid_size=redshift_grid_size,
            mass_grid_size=mass_grid_size,
            purity=purity_aguena,
            completeness=completeness_aguena
        )
        recipe_grid_abundance.setup()

        ssc_cov    = cov_terms["SSC"].copy()
        counts_cov = cov_terms["gauss"].copy()

        data_type = sacc.standard_types.cluster_counts
        data_points = sacc_full.get_data_points(data_type=data_type)
        theory_counts = {}

        for d_point in data_points:
            trs = d_point.tracers
            point_idx = sacc_full.indices(data_type=data_type, tracers=trs)
            tr_objs = [sacc_full.get_tracer(tr) for tr in trs]
            area   = tr_objs[0].sky_area
            rich_l = tr_objs[1].lower
            rich_u = tr_objs[1].upper
            z_l    = tr_objs[2].lower
            z_u    = tr_objs[2].upper
            k = point_idx[0]
            zb = [z_l, z_u]
            pb = [rich_l, rich_u]
            theory_counts[k] = recipe_grid_abundance.evaluate_theory_prediction_counts(
                np.array(zb), np.array(pb), area
            )

        for i in theory_counts:
            tjpcov_counts_i = counts_cov[i, i]
            ssc_term_i      = ssc_cov[i, i]
            old_cov_ii = full_cov[i, i]
            full_cov[i, i] = theory_counts[i] + ssc_term_i * theory_counts[i]**2 / (tjpcov_counts_i**2)
            print(f"Replaced cov points at {i, i}. From {old_cov_ii} to {full_cov[i, i]}")
            for j in theory_counts:
                if j <= i:
                    continue
                old_cov_ij      = full_cov[i, j]
                tjpcov_counts_j = counts_cov[j, j]
                ssc_term_ij     = ssc_cov[i, j]
                val = ssc_term_ij * theory_counts[i] * theory_counts[j] / (tjpcov_counts_i * tjpcov_counts_j)
                full_cov[i, j] = val
                full_cov[j, i] = val
                print(f"Replaced cov points at {i, j}. From {old_cov_ij} to {full_cov[i, j]}")

        return full_cov