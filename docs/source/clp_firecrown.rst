==================================================
CLPFirecrown Configuration Options
==================================================

This document describes the configuration options available for the
CLPFirecrown stage. These options control the construction of the
Firecrown likelihood, the cluster modeling assumptions, and the sampling
configuration.

General Options
===============

``hmf`` (str, default: "bocquet16")
    Halo Mass Function used in the analysis.
    Supported values:

    * angulo12
    * bocquet16
    * bocquet20
    * despali16
    * jenkins01
    * press74
    * sheth99
    * tinker08
    * tinker10
    * watson13

        m_pivot (float)
        z_pivot (float)

    Richness–mass relation parameters:
        mu_p0, mu_p1, mu_p2
        sigma_p0, sigma_p1, sigma_p2

    These control the mapping between halo mass and observable richness.

--------------------------------------------------
Pipeline-Specific Options
--------------------------------------------------

replace_tjpcov_cov (bool, default: True)
    If True, replaces the TJPCov cluster-count covariance with a
    custom covariance computed using CROW.

    This includes:
        - Recomputing diagonal terms using theory predictions
        - Combining SSC and Gaussian contributions

    WARNING:
        This is a temporary workaround and should ideally be handled
        directly inside TJPCov.

wazp_catalog (bool, default: False)
    Special handling for WaZP catalogs:
        - Disables purity model
        - Modifies completeness parameters

--------------------------------------------------
Internal Behavior
--------------------------------------------------

1. Covariance Computation
-------------------------

TJPCov computes covariance terms based on cov_type and writes:

    clusters_sacc_file_cov.sacc

Additional intermediate files may include:
    - clusters_sacc_file_cov_SSC.sacc
    - clusters_sacc_file_cov_gauss.sacc

2. Mixed Data Handling
----------------------

If the input SACC file contains multiple observables:

    - Cluster counts covariance is recomputed
    - Other covariance blocks (e.g. lensing) are preserved

This is handled by:
    extract_data_covariance()

3. Custom Covariance Replacement
--------------------------------

If replace_tjpcov_cov = True:

    - The covariance is recomputed using CROW
    - Theory predictions for counts are evaluated
    - Covariance elements are replaced using:

        C_ii = N_theory + SSC * (N_theory^2 / N_gauss^2)
        C_ij = SSC_ij * (N_i * N_j) / (N_gauss_i * N_gauss_j)

    This modifies only the cluster-count block.

--------------------------------------------------
Known Limitations
--------------------------------------------------

- replace_crow_counts() is a temporary workaround and should not exist
  long-term (should be implemented in TJPCov directly).

- Some configuration options (e.g. cosmology) are duplicated and not
  consistently used across all steps.

- Hardcoded choices exist:
    - Mass function (Despali16)
    - Grid sizes
    - Recipe type (GridBinnedClusterRecipe)

- No validation is performed on input configuration.

--------------------------------------------------
Example Configuration
--------------------------------------------------

CLPCovariance:
    use_mpi: False
    do_xi: False
    cov_type: [ClusterCountsGaussian, ClusterCountsSSC]

    cosmo: 'set'

    parameters:
        Omega_c: 0.22
        Omega_b: 0.0448
        h: 0.71
        n_s: 0.963
        sigma8: 0.8
        w0: -1
        wa: 0
        transfer_function: 'boltzmann_camb'

    photo-z:
        sigma_0: 0.05

    mor_parameters:
        mass_func: 'Despali16'
        mass_def: '200c'
        halo_bias: 'Tinker10'
        min_halo_mass: 1.0e12
        max_halo_mass: 3.16e15

        m_pivot: 14.3
        z_pivot: 0.5

        mu_p0: 3.3439
        mu_p1: 0.9582
        mu_p2: -0.0193
        sigma_p0: 0.5623
        sigma_p1: 0.0455
        sigma_p2: -0.0445

--------------------------------------------------

Notes:
- Additional TJPCov parameters can be added depending on the covariance model.
- Ensure consistency between SACC input data and MOR configuration.

Firecrown Likelihood Parameters
===============================

``firecrown_parameters`` (dict)
    Parameters controlling the cluster mass–observable relation and systematics.

    Same format as cosmological parameters:

    - sample: True/False
    - values: [min, fiducial, max] or scalar

    Typical parameters include:

    - mass_distribution_mu*
    - mass_distribution_sigma*
    - completeness_*
    - purity_*
    - cluster_theory_cluster_concentration

    Additional parameters can be added freely as long as they are
    recognized by the Firecrown likelihood.

Notes and Caveats
=================

- If both ``use_cluster_counts`` and ``use_shear_profile`` are True,
  the likelihood combines both observables.

- If only one is enabled, only that observable is used.

- Purity can be forcibly disabled internally when using shear-only
  configurations.

- The SACC file provided as input must match the ``survey_name``.

- Grid-based recipes are recommended for performance unless exact
  integration is required.

- Some parameters (e.g. ``cluster_concentration``) are passed directly
  to the underlying modeling code and are not validated here.

Example Configuration
=====================

.. code-block:: yaml

    CLPFirecrown:
        hmf: 'despali16'
        min_mass: 12.0
        max_mass: 15.5
        min_z: 0.2
        max_z: 0.8
        mass_def: '200c'
        use_shear_profile: True
        use_completeness: True
        use_purity: True
        use_grid: True
        is_deltasigma: True
        use_beta_interp: False
        beta_parameters: [10.0, 5.0]
        pivot_mass: 14.3
        pivot_z: 0.5
        survey_name: 'cosmodc2_redmapper'
        sampler: 'emcee'
        use_cluster_counts: True
        use_mean_log_mass: False
        use_mean_deltasigma: True
        emcee_walkers: 50
        emcee_samples: 20000
        emcee_nsteps: 20

        

Firecrown parameters can be extended as needed depending on the
likelihood model.