==================================================
CLPFirecrown Configuration Options
==================================================

This document describes the configuration options of the CLPFirecrown
stage. The stage prepares a cluster cosmology inference with Firecrown and
CosmoSIS: it writes the Firecrown likelihood and the CosmoSIS configuration
files. It does not run the sampler. Sampling is a separate CosmoSIS run
(see `Running the Inference`_).

Overview
========

CLPFirecrown runs the following steps:

1. Read the SACC file with covariance (from CLPCovariance) and the fiducial
   cosmology.
2. Write the Firecrown likelihood file. It builds CROW recipes for the
   cluster counts and/or the stacked shear profile, from the modeling
   options below.
3. Write the CosmoSIS pipeline file: consistency, CAMB and the Firecrown
   likelihood modules, plus the sampler settings.
4. Write the CosmoSIS values file: cosmological and Firecrown parameters,
   with their priors.

Inputs and outputs:

- ``clusters_sacc_file_cov`` (input): SACC data vector with covariance, from
  CLPCovariance
- ``fiducial_cosmology`` (input): fiducial cosmology, shared with TXPipe and
  CLPCovariance
- ``sampler_file`` (output): CosmoSIS pipeline file (``sampler_file.ini``)
- ``likelihood_file`` (output): Firecrown likelihood (``likelihood_file.py``)
- ``priors_file`` (output): CosmoSIS values file with priors
  (``priors_file.ini``)

Theory Model
============

The predictions are computed with CROW. For a redshift bin :math:`i` and a
richness bin :math:`j`, the number counts and the stacked lensing profile
are

.. math::

    N_{ij} = {} & \Omega_S \int_{z_i}^{z_{i+1}} dz \int_{\ln\lambda_j}^{\ln\lambda_{j+1}} d\ln\lambda \int_{M_{\min}}^{M_{\max}} dM \\
    & \times \frac{d^2V}{dz \, d\Omega} \, \frac{dn}{dM}(M, z) \, P(\ln\lambda \mid M, z) \, \Phi(M, \lambda, z)

.. math::

    \Theta_{ij}(R) = {} & \frac{\Omega_S}{N_{ij}} \int_{z_i}^{z_{i+1}} dz \int_{\ln\lambda_j}^{\ln\lambda_{j+1}} d\ln\lambda \int_{M_{\min}}^{M_{\max}} dM \\
    & \times \frac{d^2V}{dz \, d\Omega} \, \frac{dn}{dM}(M, z) \, P(\ln\lambda \mid M, z) \, \Phi(M, \lambda, z) \, \Theta(R \mid M)

where:

- :math:`\Omega_S`: survey area in steradians, from the SACC file
- :math:`dn/dM`: halo mass function (``hmf``, ``mass_def``)
- :math:`P(\ln\lambda \mid M, z)`: mass–richness relation (see below)
- :math:`\Phi = c(M, z) / p(\lambda, z)`: selection function, from the
  completeness :math:`c` and the purity :math:`p` (see below). Each is 1
  when disabled.
- :math:`\Theta(R \mid M)`: lensing profile of a halo of mass :math:`M`,
  either the excess surface density :math:`\Delta\Sigma` or the reduced
  tangential shear :math:`g_t`

Mass–richness relation
----------------------

The richness follows the
`Murata et al. (2019) <https://doi.org/10.1093/pasj/psz092>`__ model. At
fixed mass and redshift, :math:`\ln\lambda` is Gaussian. Its mean and its
scatter are each linear in :math:`\ln M` and :math:`\ln(1+z)`, with 3
parameters each:

.. math::

    P(\ln\lambda \mid M, z) = \frac{1}{\sqrt{2\pi} \, \sigma_{\ln\lambda}} \exp\left[ -\frac{(\ln\lambda - \mu_{\ln\lambda})^2}{2 \sigma_{\ln\lambda}^2} \right]

.. math::

    \mu_{\ln\lambda}(M, z) = \mu_0 + \mu_m \ln\frac{M}{M_{\rm piv}} + \mu_z \ln\frac{1+z}{1+z_{\rm piv}}

.. math::

    \sigma_{\ln\lambda}(M, z) = \sigma_0 + \sigma_m \ln\frac{M}{M_{\rm piv}} + \sigma_z \ln\frac{1+z}{1+z_{\rm piv}}

The pivots are set by ``pivot_mass`` and ``pivot_z``. The 6 parameters are
Firecrown parameters (see `Firecrown Parameters`_).

Completeness
------------

With ``use_completeness``, the
`Aguena & Lima (2018) <https://arxiv.org/abs/1611.05468>`__ completeness
model:

.. math::

    c(M, z) = \frac{\left(M / M_0(z)\right)^{n_c(z)}}{1 + \left(M / M_0(z)\right)^{n_c(z)}}

.. math::

    n_c(z) = a_n + b_n \, (1 + z), \qquad \log_{10} M_0(z) = a_{\rm piv} + b_{\rm piv} \, (1 + z)

Purity
------

With ``use_purity``, the Aguena & Lima (2018) purity model, written in
:math:`\ln\lambda`:

.. math::

    p(\lambda, z) = \frac{\left(\ln\lambda / \ln\lambda_0(z)\right)^{n_p(z)}}{1 + \left(\ln\lambda / \ln\lambda_0(z)\right)^{n_p(z)}}

.. math::

    n_p(z) = a_n + b_n \, (1 + z), \qquad \ln\lambda_0(z) = a_{\rm piv} + b_{\rm piv} \, (1 + z)

Lensing profile
---------------

:math:`\Theta(R \mid M)` is computed with CLMM from an NFW profile (one-halo
term). The concentration is the Firecrown parameter
``cluster_theory_cluster_concentration``. A value below 1 uses the
`Bhattacharya et al. (2013) <https://arxiv.org/abs/1112.5479>`__
mass–concentration relation instead. The two-halo term and the boost-factor
correction can be added (``two_halo_term``, ``boost_factor``).

For the reduced shear :math:`g_t`, the lensing efficiency is averaged over
the source redshift distribution (``beta_parameters``, ``use_beta_interp``).

Likelihood
----------

The likelihood is a Gaussian (Firecrown ``ConstGaussian``) over all the
selected statistics:

- ``BinnedClusterNumberCounts``: counts :math:`N_{ij}` (and optionally the
  mean log mass)
- ``BinnedClusterShearProfile``: stacked :math:`\Delta\Sigma` or :math:`g_t`
  profiles

The covariance is read from the SACC file and held fixed during sampling.

Modeling Options
================

``hmf`` (str, default: "despali16")
    Halo mass function (case-insensitive). Supported values:

    - angulo12
    - bocquet16
    - bocquet20
    - despali16
    - jenkins01
    - press74
    - sheth99
    - tinker08
    - tinker10
    - watson13

``mass_def`` (str, default: "200c")
    Mass definition used by the halo mass function (e.g. ``"200c"``,
    ``"500c"``).

``min_mass``, ``max_mass`` (float, default: 12.0, 15.5)
    Halo mass integration range, in :math:`\log_{10}(M / M_\odot)`.
    Note: CLPCovariance uses linear masses for the same range.

``min_z``, ``max_z`` (float, default: 0.2, 0.8)
    True redshift integration range.

``pivot_mass`` (float, default: 14.3)
    Pivot mass of the mass–richness relation,
    :math:`\log_{10}(M_{\rm piv} / M_\odot)`.

``pivot_z`` (float, default: 0.5)
    Pivot redshift of the mass–richness relation.

``survey_name`` (str, default: "cosmodc2_redmapper")
    Survey tracer name in the SACC file. Must match it.

Observable Options
==================

``use_cluster_counts`` (bool, default: True)
    Include the cluster number counts in the likelihood.

``use_shear_profile`` (bool, default: False)
    Include the stacked shear profiles in the likelihood.

``is_deltasigma`` (bool, default: False)
    Theory side: predict :math:`\Delta\Sigma` (True) or the reduced shear
    :math:`g_t` (False).

``use_mean_deltasigma`` (bool, default: False)
    Data side: read :math:`\Delta\Sigma` (True) or the reduced shear
    :math:`g_t` (False) from the SACC file. Set it to the same value as
    ``is_deltasigma``.

``use_mean_log_mass`` (bool, default: False)
    Also use the mean log mass per bin, read from the SACC file.

Selection Function Options
==========================

``use_completeness`` (bool, default: True)
    Apply the completeness model :math:`c(M, z)` (see `Completeness`_).

``use_purity`` (bool, default: True)
    Apply the purity model :math:`p(\lambda, z)` (see `Purity`_).

Lensing Options
===============

Only used when ``use_shear_profile`` is True.

``two_halo_term`` (bool, default: False)
    Add the two-halo term to the lensing profile.

``boost_factor`` (bool, default: False)
    Apply the boost-factor correction to the lensing profile.

``set_concentration`` (bool)
    Reserved for future use. The stage does not read it yet, so it
    currently has no effect. The concentration is set with
    ``cluster_theory_cluster_concentration`` in ``firecrown_parameters``.

``use_beta_interp`` (bool, default: False)
    Interpolate the mean lensing efficiency :math:`\langle\beta_s\rangle` in
    redshift. Only used for the reduced shear (``is_deltasigma: False``).

``beta_parameters`` (list of float, default: [10.0, 5.0])
    Parameters of the mean lensing efficiency :math:`\langle\beta_s\rangle`,
    passed to CROW as ``[z_inf, zmax]``: the redshift taken as infinity, and
    the maximum source redshift. Only used for the reduced shear
    (``is_deltasigma: False``).

    The first value (``z_inf``, 10 by default) is also always used as the
    maximum redshift of CAMB (``zmax`` in the ``[camb]`` section of
    ``sampler_file.ini``), even without shear profiles.

Integration Options
===================

``use_grid`` (bool, default: True)
    If True, use ``GridBinnedClusterRecipe``: integrals on precomputed
    grids, about 100 times faster, recommended for inference.
    If False, use ``ExactBinnedClusterRecipe`` (direct integration with
    NumCosmo), slower, useful to validate the grid results.

``redshift_grid_size`` (int, default: 20)
    Number of redshift grid points.

``mass_grid_size`` (int, default: 60)
    Number of mass grid points.

``proxy_grid_size`` (int, default: 20)
    Number of richness grid points.

The grid sizes are only used when ``use_grid`` is True.

Sampler Options
===============

``sampler`` (str, default: "emcee")
    CosmoSIS sampler. Supported values:

    - ``test``: a single likelihood evaluation, useful to check the setup
    - ``metropolis``
    - ``emcee``
    - ``polychord``

``emcee_walkers`` (int, default: 100)
    Number of emcee walkers.

``emcee_samples`` (int, default: 20000)
    Number of emcee samples per walker.

``emcee_nsteps`` (int, default: 20)
    Number of emcee steps between chain outputs.

``polychord_live_points`` (int, default: 500)
    Number of PolyChord live points.

``polychord_num_repeats`` (int, default: 30)
    Number of PolyChord slice-sampling repeats.

``polychord_tolerance`` (float, default: 0.05)
    PolyChord evidence tolerance (stopping criterion).

``polychord_feedback`` (int, default: 1)
    PolyChord verbosity level.

.. note::

    The ``polychord_*`` options fill the ``[polychord]`` section of
    ``sampler_file.ini`` (see the
    `CosmoSIS PolyChord sampler <https://cosmosis.readthedocs.io/en/latest/reference/samplers/polychord.html>`__).

``resume`` (bool, default: False)
    If True, CosmoSIS appends to the existing chain instead of starting a
    new one.

``filename`` (str, optional, default: "output_rp/number_counts_samples.txt")
    Path of the output chain file.

Cosmological Parameters
=======================

``tau`` (float, default: 0.08)
    Optical depth to reionization. It is a CAMB input, not part of the
    fiducial cosmology file, so it is set here.

``cosmological_parameters`` (dict, default: {})
    Overrides on top of the fiducial cosmology. By default, every
    cosmological parameter is fixed to its value in ``fiducial_cosmology``.
    List a parameter here to sample it.

    Parameter names (CosmoSIS naming): ``omega_c``, ``omega_b``, ``h0``,
    ``n_s``, ``sigma_8``, ``omega_k``, ``w``, ``wa``, ``tau``.

    A fixed entry (``sample: False``) must have the fiducial value, otherwise
    the stage fails. To change the fiducial cosmology, edit the
    ``fiducial_cosmology`` file instead, so that all stages use it.

    Example:

    .. code-block:: yaml

        cosmological_parameters:
            omega_c: {'sample': True, 'values': [0.10, 0.22, 0.5]}
            sigma_8: {'sample': True, 'values': [0.5, 0.800, 1.1]}

Firecrown Parameters
====================

``firecrown_parameters`` (dict, default: {})
    Parameters of the CROW models, written to the
    ``[firecrown_number_counts]`` section of the values file. Set all the
    parameters of the models you enable, fixed or sampled.

    **Mass–richness relation:**

    - ``mass_distribution_mu0``, ``mass_distribution_mu1``,
      ``mass_distribution_mu2``: :math:`\mu_0, \mu_m, \mu_z`
    - ``mass_distribution_sigma0``, ``mass_distribution_sigma1``,
      ``mass_distribution_sigma2``: :math:`\sigma_0, \sigma_m, \sigma_z`

    **Completeness** (``use_completeness``):

    - ``completeness_a_n``, ``completeness_b_n``: :math:`a_n, b_n`
    - ``completeness_a_logm_piv``, ``completeness_b_logm_piv``:
      :math:`a_{\rm piv}, b_{\rm piv}`
    - CROW defaults (cosmoDC2 redMaPPer): 1.1321, 0.7751, 13.31, 0.2025

    **Purity** (``use_purity``):

    - ``purity_a_n``, ``purity_b_n``: :math:`a_n, b_n`
    - ``purity_a_logm_piv``, ``purity_b_logm_piv``:
      :math:`a_{\rm piv}, b_{\rm piv}`
    - CROW defaults (cosmoDC2 redMaPPer): 1.9830, 0.8121, 2.2183, -0.6592

    **Lensing profile** (``use_shear_profile``):

    - ``cluster_theory_cluster_concentration``: halo concentration
      (below 1: Bhattacharya et al. 2013 relation)

Parameter Format
================

Entries in ``cosmological_parameters`` and ``firecrown_parameters`` use the
same format:

.. code-block:: yaml

    # Sampled, with a flat prior between min and max:
    name: {'sample': True, 'values': [min, start, max]}

    # Fixed:
    name: {'sample': False, 'values': value}

Pipeline Configuration
======================

The stage is wired into a ceci pipeline file. This example comes from
``examples/cosmodc2_redmapper/baseline/cosmodc2_redmapper_full_analysis/run_in2p3_both/Firecrown.yml``:

.. code-block:: yaml

    id: Firecrown
    modules: clpipe
    launcher:
        name: mini
        interval: 0.5
    site:
        name: local
        max_threads: 4
    stages:
      - name: CLPFirecrown
        module_name: clpipe.clp_firecrown
        nprocess: 1
    inputs:
        fiducial_cosmology: /sps/lsst/groups/clusters/cl_pipeline_project/TXPipe_data/cosmodc2/fiducial_cosmology.yml
        clusters_sacc_file_cov: ./outputs_both/clusters_sacc_file_cov.sacc
    config: ./config_in2p3_both.yml
    resume: false
    output_dir: ./outputs_both
    log_dir: ./logs_both

The stage options go in the stage config file (``config:`` above), under a
``CLPFirecrown`` block.

Running the Inference
=====================

Generate the files with ceci, then run CosmoSIS from the output directory:

.. code-block:: bash

    ceci Firecrown.yml

    cd ./outputs_both
    cosmosis sampler_file.ini

    # or, with MPI:
    mpirun -n 30 cosmosis --mpi sampler_file.ini

The generated files refer to each other, and to the SACC file, by file name
only. Run CosmoSIS from the output directory, with the SACC file in it (this
is the case when CLPCovariance writes to the same output directory).
CosmoSIS also needs the ``CSL_DIR`` environment variable (see the
installation instructions in the README).

Known Limitations
=================

- The ``test`` and ``metropolis`` sampler settings are fixed in the
  generated file (``metropolis``: 1000 samples).

- The covariance is held fixed during sampling.

- Only the binned likelihood is supported.

- The names in ``firecrown_parameters`` are not validated by the stage.

Example Configuration
=====================

This is the ``CLPFirecrown`` block of the baseline cosmoDC2 redMaPPer
analysis (counts + stacked :math:`\Delta\Sigma`):

.. code-block:: yaml

    CLPFirecrown:
        hmf: 'despali16'
        mass_def: '200c'
        min_mass: 12.0
        max_mass: 15.5
        min_z: 0.2
        max_z: 0.8
        pivot_mass: 14.3
        pivot_z: 0.5
        survey_name: 'cosmodc2_redmapper'

        use_cluster_counts: true
        use_shear_profile: true
        is_deltasigma: true
        use_mean_deltasigma: true
        use_mean_log_mass: false

        use_completeness: true
        use_purity: false

        use_grid: true
        use_beta_interp: false
        beta_parameters: [10.0, 5.0]

        sampler: 'emcee'
        emcee_walkers: 100
        emcee_samples: 50000
        emcee_nsteps: 20

        cosmological_parameters:
            omega_c: {'sample': true, 'values': [0.10, 0.22, 0.5]}
            sigma_8: {'sample': true, 'values': [0.5, 0.800, 1.1]}

        firecrown_parameters:
            mass_distribution_mu0: {'sample': true, 'values': [2.0, 3.3439, 10.0]}
            mass_distribution_mu1: {'sample': true, 'values': [0.5, 0.958236982, 2.0]}
            mass_distribution_mu2: {'sample': true, 'values': [-2.0, -0.0192802, 2.0]}
            mass_distribution_sigma0: {'sample': true, 'values': [0.1, 0.562317194, 2.0]}
            mass_distribution_sigma1: {'sample': true, 'values': [-0.6, 0.04552506, 0.3]}
            mass_distribution_sigma2: {'sample': true, 'values': [-0.5, -0.0445, 2.0]}
            completeness_a_n: {'sample': false, 'values': 1.1321}
            completeness_b_n: {'sample': false, 'values': 0.7751}
            completeness_a_logm_piv: {'sample': false, 'values': 13.31}
            completeness_b_logm_piv: {'sample': false, 'values': 0.2025}
            cluster_theory_cluster_concentration: {'sample': false, 'values': 3.8}

References
==========

- Murata et al. (2019), mass–richness relation:
  `doi:10.1093/pasj/psz092 <https://doi.org/10.1093/pasj/psz092>`__
  (`arXiv:1904.07524 <https://arxiv.org/abs/1904.07524>`__)
- Aguena & Lima (2018), completeness and purity:
  `arXiv:1611.05468 <https://arxiv.org/abs/1611.05468>`__
- Despali et al. (2016), halo mass function:
  `arXiv:1507.05627 <https://arxiv.org/abs/1507.05627>`__
- Navarro, Frenk & White (1997), NFW profile:
  `arXiv:astro-ph/9611107 <https://arxiv.org/abs/astro-ph/9611107>`__
- Bhattacharya et al. (2013), mass–concentration relation:
  `arXiv:1112.5479 <https://arxiv.org/abs/1112.5479>`__
- Zuntz et al. (2015), CosmoSIS:
  `arXiv:1409.3409 <https://arxiv.org/abs/1409.3409>`__
- Foreman-Mackey et al. (2013), emcee:
  `arXiv:1202.3665 <https://arxiv.org/abs/1202.3665>`__
- Firecrown: https://github.com/LSSTDESC/firecrown
- CROW: https://github.com/LSSTDESC/crow
