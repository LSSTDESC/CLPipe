#!/usr/bin/env python
#File dedicated to impement the Firecrown pipeline stage into ceci
from .file_types import SACCFile, PythonFile, CosmosisFile, FiducialCosmology
from ceci import PipelineStage
from ceci.config import StageParameter

import sys
import os
import shutil
import re
import textwrap
import math
import numpy as np

# CosmoSIS's [cosmological_parameters] naming differs from the fiducial
# cosmology file's CCL-native naming. Maps CCL key -> CosmoSIS key.
_CCL_TO_COSMOSIS_COSMO_MAP = {
    "Omega_c": "omega_c",
    "Omega_b": "omega_b",
    "h": "h0",
    "n_s": "n_s",
    "sigma8": "sigma_8",
    "Omega_k": "omega_k",
    "w0": "w",
    "wa": "wa",
}


class CLPFirecrown(PipelineStage):
    """Firecrown pipeline stage for cluster cosmology analysis.

    This stage:
    - Builds a Firecrown likelihood from a SACC file
    - Generates the corresponding CosmoSIS configuration
    - Writes parameter files for sampling, defaulting the cosmological
      block to the fiducial cosmology (shared with TXPipe/TJPCov); any
      entry in `cosmological_parameters` overrides that default, which is
      how a parameter gets sampled instead of held fixed

    Key configuration groups:
    - Modeling options (hmf, mass range, redshift range)
    - Observable selection (cluster counts, shear)
    - Systematics (purity, completeness)
    - Sampling configuration (emcee, polychord)

    Full configuration documentation:
    See docs/clp_firecrown.txt
    """
    name = "CLPFirecrown"

    inputs = [
        ("clusters_sacc_file_cov", SACCFile),  # For firecrown Likelihood
        ("fiducial_cosmology", FiducialCosmology),
    ]

    outputs = [
        ("sampler_file", CosmosisFile),
        ("likelihood_file", PythonFile),
        ("priors_file", CosmosisFile),
    ]

    config_options = {
        # Modeling
        "hmf": StageParameter(str, "despali16", msg="Halo mass function model."),
        "min_mass": StageParameter(float, 12.0, msg="Minimum log10 halo mass."),
        "max_mass": StageParameter(float, 15.5, msg="Maximum log10 halo mass."),
        "min_z": StageParameter(float, 0.2, msg="Minimum cluster redshift."),
        "max_z": StageParameter(float, 0.8, msg="Maximum cluster redshift."),
        "mass_def": StageParameter(str, "200c", msg="Halo mass definition."),
        "pivot_mass": StageParameter(float, 14.3, msg="Pivot log10 halo mass for the mass-richness relation."),
        "pivot_z": StageParameter(float, 0.5, msg="Pivot redshift for the mass-richness relation."),
        "survey_name": StageParameter(str, "cosmodc2_redmapper", msg="Survey name used in the SACC file."),
        # Observable selection
        "use_cluster_counts": StageParameter(bool, True, msg="Include cluster number counts."),
        "use_mean_log_mass": StageParameter(bool, False, msg="Include mean log-mass observable."),
        "use_mean_deltasigma": StageParameter(bool, False, msg="Include DeltaSigma observable."),
        "use_shear_profile": StageParameter(bool, False, msg="Generate the cluster shear profile likelihood."),
        # Systematics
        "use_completeness": StageParameter(bool, True, msg="Apply completeness model."),
        "use_purity": StageParameter(bool, True, msg="Apply purity model."),
        "use_grid": StageParameter(bool, True, msg="Use the gridded recipe."),
        "is_deltasigma": StageParameter(bool, False, msg="Use DeltaSigma instead of reduced shear."),
        "use_beta_interp": StageParameter(bool, False, msg="Use beta interpolation."),
        "beta_parameters": StageParameter(list, [10.0, 5.0], msg="Parameters for beta interpolation."),
        "redshift_grid_size": StageParameter(int, 20, msg="Number of redshift grid points."),
        "mass_grid_size": StageParameter(int, 60, msg="Number of mass grid points."),
        "proxy_grid_size": StageParameter(int, 20, msg="Number of richness grid points."),
        "two_halo_term": StageParameter(bool, False, msg="Include the two-halo term."),
        "boost_factor": StageParameter(bool, False, msg="Apply the boost-factor correction."),
        # Sampling
        "sampler": StageParameter(str, "emcee", msg="CosmoSIS sampler."),
        "emcee_walkers": StageParameter(int, 100, msg="Number of emcee walkers."),
        "emcee_samples": StageParameter(int, 20000, msg="Number of emcee samples."),
        "emcee_nsteps": StageParameter(int, 20, msg="Number of emcee steps per sample."),
        "polycord_live_points": StageParameter(int, 500, msg="Number of PolyChord live points."),
        "polycord_num_repeats": StageParameter(int, 30, msg="Number of PolyChord repeats."),
        "polycord_tolerance": StageParameter(float, 0.05, msg="PolyChord evidence tolerance."),
        "polycord_feedback": StageParameter(int, 1, msg="PolyChord feedback level."),
        "resume": StageParameter(bool, False, msg="If True, CosmoSIS appends to the existing chain at `filename` instead of starting fresh"),
        # Cosmology -- tau is a CAMB input (reionization optical depth), not
        # a CCL cosmology parameter, so it isn't in the fiducial cosmology
        # file and has to stay a stage config option.
        "tau": StageParameter(float, 0.08, msg="CMB optical depth to reionization (CAMB input, not part of the CCL fiducial cosmology)."),
        # Parameter blocks
        "cosmological_parameters": StageParameter(
            dict, {},
            msg=(
                "Overrides on top of the fiducial cosmology. Only needs to "
                "contain entries you want to sample or fix to a non-fiducial "
                "value; anything not listed here is taken from "
                "fiducial_cosmology as a fixed value."
            ),
        ),
        "firecrown_parameters": StageParameter(dict, {}, msg="Dictionary describing Firecrown likelihood parameters."),
    }

    def run(self):
        """Run the analysis for this stage.
        Generates Firecrown likelihood and
        cosmosis ini files.
        """

        output_cosmosis_file = self.get_output(
            'sampler_file',
            final_name=True
        )
        output_likelihood_file = self.get_output(
            'likelihood_file',
            final_name=True
        )
        output_parameters_file = self.get_output(
            'priors_file',
            final_name=True
        )

        ok_python = self.generate_python_file(output_likelihood_file)
        ok_ini = self.generate_ini_file(
            output_cosmosis_file,
            os.path.basename(output_likelihood_file),
            os.path.basename(output_parameters_file),
        )
        ok_params = self.generate_cosmosis_parameters_file(output_parameters_file)

        if not (ok_python and ok_ini and ok_params):
            raise RuntimeError(
                "CLPFirecrown: one or more output files failed to "
                "generate. See the printed errors above for details. "
                f"(python={ok_python}, ini={ok_ini}, params={ok_params})"
            )

    def _fiducial_cosmological_parameters(self):
        """Fiducial cosmology (shared with TXPipe), translated into
        CosmoSIS's [cosmological_parameters] naming and CosmoSIS-parameter
        dict shape, all fixed (sample=False). tau has no fiducial-file
        equivalent and comes from the stage config instead.
        """
        with self.open_input("fiducial_cosmology", wrapper=True) as f:
            raw = f.content

        params = {
            cosmosis_name: {"sample": False, "values": raw[ccl_name]}
            for ccl_name, cosmosis_name in _CCL_TO_COSMOSIS_COSMO_MAP.items()
        }
        params["tau"] = {"sample": False, "values": self.config["tau"]}
        return params

    def generate_python_file(self, path_name):
        """Generates a Python file based on the configuration dictionary.

        Args:
            path_name (str): Path to save the generated Python file.
        """
        my_configs = self.config
        hmf_dict = {
            'angulo12': 'ccl.halos.MassFuncAngulo12',
            'bocquet16': 'ccl.halos.MassFuncBocquet16',
            'bocquet20': 'ccl.halos.MassFuncBocquet20',
            'despali16': 'ccl.halos.MassFuncDespali16',
            'jenkins01': 'ccl.halos.MassFuncJenkins01',
            'press74': 'ccl.halos.MassFuncPress74',
            'sheth99': 'ccl.halos.MassFuncSheth99',
            'tinker08': 'ccl.halos.MassFuncTinker08',
            'tinker10': 'ccl.halos.MassFuncTinker10',
            'watson13': 'ccl.halos.MassFuncWatson13',
        }
        try:
            cfg = self.config

            # Input file
            sacc_path = self.get_input("clusters_sacc_file_cov")
            sacc_filename = os.path.basename(sacc_path)

            # Halo mass function
            hmf_key = cfg["hmf"].lower()
            try:
                hmf = hmf_dict[hmf_key]
            except KeyError:
                raise ValueError(
                    f"Unknown halo mass function '{hmf_key}'. "
                    f"Available options are: {', '.join(sorted(hmf_dict))}"
                )

            # Modeling
            mass_def = cfg["mass_def"]
            min_mass = cfg["min_mass"]
            max_mass = cfg["max_mass"]
            min_z = cfg["min_z"]
            max_z = cfg["max_z"]
            pivot_mass = cfg["pivot_mass"]
            pivot_z = cfg["pivot_z"]
            survey_name = cfg["survey_name"]

            # Observables
            use_cluster_counts = cfg["use_cluster_counts"]
            use_shear_profile = cfg["use_shear_profile"]
            use_mean_log_mass = cfg["use_mean_log_mass"]
            use_mean_deltasigma = cfg["use_mean_deltasigma"]

            # Systematics
            use_completeness = cfg["use_completeness"]
            use_purity = cfg["use_purity"]

            # Grid / recipe options
            use_grid = cfg["use_grid"]
            is_deltasigma = cfg["is_deltasigma"]
            use_beta_interp = cfg["use_beta_interp"]
            beta_parameters = cfg["beta_parameters"]
            redshift_grid_size = cfg["redshift_grid_size"]
            mass_grid_size = cfg["mass_grid_size"]
            proxy_grid_size = cfg["proxy_grid_size"]
            two_halo_term = cfg["two_halo_term"]
            boost_factor = cfg["boost_factor"]
            # Open the file to be written
            with open(path_name, "w") as f:
                f.write("import os\n\n")

                f.write("import pyccl as ccl\n")
                f.write("import sacc\n\n")

                # Core CROW imports
                f.write("from crow import ClusterAbundance, ClusterShearProfile, kernel, mass_proxy\n")
                f.write("from crow.properties import ClusterProperty\n")
                f.write("from crow.recipes.binned_grid import GridBinnedClusterRecipe\n\n")
                f.write(
                    "from crow import purity_models, completeness_models\n"
                )
                # Firecrown likelihoods
                f.write(
                    "from firecrown.likelihood import (\n"
                    "    ConstGaussian,\n"
                    "    BinnedClusterShearProfile,\n"
                    "    BinnedClusterNumberCounts,\n"
                    "    Likelihood,\n"
                    "    NamedParameters,\n"
                    ")\n"
                )
                f.write("from firecrown.modeling_tools import ModelingTools\n\n")
                if use_cluster_counts:
                    f.write("def get_cluster_abundance() -> ClusterAbundance:\n")
                    f.write("    \"\"\"Creates and returns a ClusterAbundance object.\"\"\" \n")
                    f.write("    cluster_theory = ClusterAbundance(\n")
                    f.write(f"    halo_mass_function = {hmf}(mass_def=\"{mass_def}\"),\n")
                    f.write("    cosmo = ccl.CosmologyVanillaLCDM()\n")
                    f.write("    )\n\n")
                    f.write("    return cluster_theory\n\n")
                if use_shear_profile:
                    f.write("\n")
                    f.write("def get_cluster_shear_profile() -> ClusterShearProfile:\n")
                    f.write("    \"\"\"Creates and returns a ClusterShearProfile object.\"\"\"\n")
                    f.write("    cluster_theory = ClusterShearProfile(\n")
                    f.write("    cosmo=ccl.CosmologyVanillaLCDM(),\n")
                    f.write(f"    halo_mass_function = {hmf}(mass_def=\"{mass_def}\"),\n")
                    f.write("    cluster_concentration=None,\n")
                    f.write(f"    is_delta_sigma={is_deltasigma},\n")
                    f.write(f"    use_beta_s_interp={use_beta_interp},\n")
                    f.write(f"    two_halo_term={two_halo_term},\n")
                    f.write(f"    boost_factor={boost_factor},\n")
                    f.write("    )\n\n")
                    f.write("    return cluster_theory\n\n")
                f.write("def get_cluster_recipe(\n")
                f.write("    cluster_theory,\n")
                f.write(f"    pivot_mass: float = {pivot_mass},\n")
                f.write(f"    pivot_redshift: float = {pivot_z},\n")
                f.write(f"    mass_interval=({min_mass}, {max_mass}),\n")
                f.write(f"    true_z_interval=({min_z}, {max_z}),\n")
                f.write(f"    is_reduced_shear = False,\n")
                f.write("):\n")
                f.write("    \"\"\"Creates and returns a ClusterRecipe.\n\n")
                f.write("    Parameters\n")
                f.write("    ----------\n")
                f.write("    cluster_theory : ClusterShearProfile or ClusterAbundance\n")
                f.write("    \"\"\"\n")
                f.write("    redshift_distribution = kernel.SpectroscopicRedshift()\n")
                if use_completeness:
                    f.write(f"    completeness = completeness_models.CompletenessAguena16()\n")
                else:
                    f.write(f"    completeness = None\n")
                if use_purity:
                    f.write(f"    purity = purity_models.PurityAguena16LnProxy()\n")
                else:
                    f.write(f"    purity = None\n")
                f.write("    if is_reduced_shear:\n")
                f.write(f"        cluster_theory.set_beta_parameters({beta_parameters[0]}, {beta_parameters[1]})\n")
                if use_beta_interp:
                    f.write(f"        cluster_theory.set_beta_s_interp(true_z_interval[0], true_z_interval[1])\n")
                if use_grid:
                    f.write(
                        "    mass_distribution = mass_proxy.MurataUnbinned(\n"
                        f"        pivot_log_mass={pivot_mass},\n"
                        f"        pivot_redshift={pivot_z},\n"
                    )
                    f.write("    )\n\n")
                    f.write("    recipe = GridBinnedClusterRecipe(\n")
                    f.write(f"        redshift_grid_size = {redshift_grid_size},\n")
                    f.write(f"        mass_grid_size = {mass_grid_size},\n")
                    f.write(f"        proxy_grid_size = {proxy_grid_size},\n")
                else:
                    f.write(
                        "    mass_distribution = mass_proxy.MurataBinned(\n"
                        f"        pivot_log_mass={pivot_mass},\n"
                        f"        pivot_redshift={pivot_z},\n"
                        "    )\n\n"
                    )
                    f.write("    recipe = ExactBinnedClusterRecipe(\n")
                f.write("        cluster_theory=cluster_theory,\n")
                f.write("        redshift_distribution=redshift_distribution,\n")
                f.write("        mass_distribution=mass_distribution,\n")
                f.write(f"        completeness=completeness,\n")
                f.write(f"        purity=purity,\n")
                f.write(f"        mass_interval=({min_mass}, {max_mass}),\n")
                f.write(f"        true_z_interval=({min_z}, {max_z}),\n")
                f.write("    )\n\n")

                f.write("    return recipe\n\n")
                f.write("def build_likelihood(build_parameters: NamedParameters) -> tuple[Likelihood, ModelingTools]:\n")
                f.write("    '''Builds the likelihood for Firecrown.''' \n")
                f.write("    # Pull params for the likelihood from build_parameters\n")
                f.write("    average_on = ClusterProperty.NONE\n")
                f.write("    if build_parameters.get_bool('use_cluster_counts', True):\n")
                f.write("        average_on |= ClusterProperty.COUNTS\n")
                f.write("    if build_parameters.get_bool('use_mean_log_mass', True):\n")
                f.write("        average_on |= ClusterProperty.MASS\n")
                f.write("    if build_parameters.get_bool('use_mean_deltasigma', True):\n")
                f.write("        average_on |= ClusterProperty.DELTASIGMA\n")
                f.write("    if build_parameters.get_bool('use_mean_reduced_shear', True):\n")
                f.write("        average_on |= ClusterProperty.SHEAR\n\n")
                f.write(f"    survey_name = '{survey_name}'\n")
                if use_shear_profile and use_cluster_counts:
                    f.write("    recipe_counts = get_cluster_recipe(get_cluster_abundance())\n")
                    f.write(f"    recipe_shear = get_cluster_recipe(get_cluster_shear_profile(), is_reduced_shear = {not is_deltasigma})\n")
                    f.write("    likelihood = ConstGaussian(\n")
                    f.write("        [\n")
                    f.write("            BinnedClusterNumberCounts(\n")
                    f.write("                average_on, survey_name, recipe_counts\n")
                    f.write("            ),\n")
                    f.write("            BinnedClusterShearProfile(\n")
                    f.write("                average_on, survey_name, recipe_shear\n")
                    f.write("            ),\n")
                    f.write("        ]\n")
                elif use_cluster_counts:
                    f.write("    recipe_counts = get_cluster_recipe(get_cluster_abundance())\n")
                    f.write("    likelihood = ConstGaussian(\n")
                    f.write("        [BinnedClusterNumberCounts(average_on, survey_name, recipe_counts)]\n")
                elif use_shear_profile:
                    f.write(f"    recipe_shear = get_cluster_recipe(get_cluster_shear_profile(), is_reduced_shear = {not is_deltasigma})\n")
                    f.write("    likelihood = ConstGaussian(\n")
                    f.write("        [BinnedClusterShearProfile(average_on, survey_name, recipe_shear)]\n")
                f.write("    )\n\n")
                f.write(f"    sacc_path = '{sacc_filename}'\n")
                f.write("    sacc_data = sacc.Sacc.load_fits(sacc_path)\n")
                f.write("    likelihood.read(sacc_data)\n\n")
                f.write("    modeling_tools = ModelingTools()\n\n")
                f.write("    return likelihood, modeling_tools\n")

            print(f"Python file generated at {path_name}")
            return True

        except Exception as e:
            print(f"Error generating file: {e}")
            return False

    def generate_ini_file(self, output_ini_path, likelihood_source_name, values_file_name):
        """Generates an .ini file.

        Args:
            output_ini_path (str): Path where the generated .ini file will be saved.
            likelihood_source_name (str): Basename of the generated Firecrown
                likelihood Python file (the `likelihood_file` output), written
                alongside output_ini_path.
            values_file_name (str): Basename of the generated CosmoSIS values
                file (the `priors_file` output), written alongside
                output_ini_path.
        """
        import cosmosis
        import firecrown

        try:
            cfg = self.config
            out_filename = cfg.get('filename', 'output_rp/number_counts_samples.txt')

            root = os.getcwd()
            FIRECROWN_DIR = os.path.dirname(firecrown.__file__)

            sampler = cfg["sampler"]
            resume = cfg["resume"]
            use_cluster_counts = cfg["use_cluster_counts"]
            use_mean_log_mass = cfg["use_mean_log_mass"]
            use_mean_deltasigma = cfg["use_mean_deltasigma"]

            emcee_walkers = cfg["emcee_walkers"]
            emcee_samples = cfg["emcee_samples"]
            emcee_nsteps = cfg["emcee_nsteps"]

            polycord_live_points = cfg["polycord_live_points"]
            polycord_num_repeats = cfg["polycord_num_repeats"]
            polycord_tolerance = cfg["polycord_tolerance"]
            polycord_feedback = cfg["polycord_feedback"]

            beta_parameters = cfg["beta_parameters"]
            with open(output_ini_path, 'w') as f:
                f.write("[runtime]\n")
                f.write(f"sampler = {sampler}\n")
                f.write(f"root = {root}\n")
                f.write(f"resume = {'T' if resume else 'F'}\n\n")
                
                f.write("[default]\n")
                f.write("fatal_errors = F\n\n")

                f.write("[output]\n")
                f.write(f"filename = {out_filename}\n")
                f.write("format = text\n")
                f.write("verbosity = 0\n\n")

                f.write("[pipeline]\n")
                f.write("modules = consistency camb firecrown_likelihood\n")
                f.write(f"values = {values_file_name}\n")
                f.write("likelihoods = firecrown\n")
                f.write("quiet = F\n")
                f.write("debug = F\n")
                f.write("timing = T\n\n")

                f.write("[consistency]\n")
                f.write("file = ${CSL_DIR}/utility/consistency/consistency_interface.py\n\n")

                f.write("[camb]\n")
                f.write("file = ${CSL_DIR}/boltzmann/camb/camb_interface.py\n\n")
                f.write("mode = all\n")
                f.write("lmax = 2500\n")
                f.write("feedback = 0\n")
                f.write("zmin = 0.0\n")
                f.write(f"zmax = {beta_parameters[0]}\n")
                f.write("nz = 100\n")
                f.write("kmin = 1e-4\n")
                f.write("kmax = 50.0\n")
                f.write("nk = 1000\n\n")

                f.write("[firecrown_likelihood]\n")
                f.write(";; Fix this to use an environment variable to find the files.\n")
                f.write(";; Set FIRECROWN_DIR to the base of the firecrown installation (or build, if you haven't installed it)\n")
                f.write(f"file = {FIRECROWN_DIR}/connector/cosmosis/likelihood.py\n")
                f.write(f"likelihood_source = {likelihood_source_name}\n")
                f.write("sampling_parameters_sections = firecrown_number_counts\n")
                f.write(f"use_cluster_counts = {str(use_cluster_counts).upper()}\n")
                f.write(f"use_mean_deltasigma = {str(use_mean_deltasigma).upper()}\n")
                f.write(f"use_mean_reduced_shear = {str(not use_mean_deltasigma).upper()}\n")
                f.write(f"use_mean_log_mass = {str(use_mean_log_mass).upper()}\n\n")
                f.write("[test]\n")
                f.write("fatal_errors = T\n")
                f.write("save_dir = output_counts_mean_mass\n\n")

                f.write("[metropolis]\n")
                f.write("samples = 1000\n")
                f.write("nsteps = 1\n\n")

                f.write("[emcee]\n")
                f.write(f"walkers = {emcee_walkers}\n")
                f.write(f"samples = {emcee_samples}\n")
                f.write(f"nsteps = {emcee_nsteps}\n")

                f.write("[polychord]\n")
                f.write(f"live_points = {polycord_live_points}\n")
                f.write(f"num_repeats = {polycord_num_repeats}\n")
                f.write(f"tolerance = {polycord_tolerance}\n")
                f.write(f"feedback = {polycord_feedback}\n")
            print(f"INI file written to {output_ini_path}")
            return True

        except Exception as e:
            print(f"Error generating INI file: {e}")
            return False

    def generate_cosmosis_parameters_file(self, output_ini_path):
        try:
            cfg = self.config
            cosmological_parameters = self._fiducial_cosmological_parameters()
            for name, override in cfg["cosmological_parameters"].items():
                if name in cosmological_parameters and not override.get("sample", False):
                    fiducial_value = cosmological_parameters[name]["values"]
                    if not math.isclose(float(override["values"]), float(fiducial_value), rel_tol=1e-6):
                        raise ValueError(
                            f"cosmological_parameters['{name}']={override['values']} "
                            f"does not match fiducial value {fiducial_value}"
                        )
                cosmological_parameters[name] = override

            with open(output_ini_path, 'w') as f:
                f.write("[cosmological_parameters]\n")
                for param, value in cosmological_parameters.items():
                    if value['sample']:
                        f.write(f"{param} = {value['values'][0]} {value['values'][1]} {value['values'][2]}\n")
                for param, value in cosmological_parameters.items():
                    if not value['sample']:
                        f.write(f"{param} = {value['values']}\n")

                f.write("[firecrown_number_counts]\n")
                for param, value in cfg["firecrown_parameters"].items():
                    if value['sample']:
                        f.write(f"{param} = {float(value['values'][0])} {float(value['values'][1])} {float(value['values'][2])}\n")
                    else:
                        f.write(f"{param} = {float(value['values'])}\n")

            print(f"Parameters INI file written to {output_ini_path}")
            return True

        except Exception as e:
            print(f"Error generating INI file: {e}")
            return False