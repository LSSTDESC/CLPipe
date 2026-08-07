from crow import ClusterShearProfile, ClusterAbundance, kernel, mass_proxy
from crow.properties import ClusterProperty
from crow.recipes.binned_grid import GridBinnedClusterRecipe
from crow.recipes.binned_exact import ExactBinnedClusterRecipe
from crow import purity_models, completeness_models
import sacc
import yaml
import pyccl as ccl
import numpy as np
from typing import Dict, Optional, Tuple, Any
import logging
from firecrown.likelihood import (
    ConstGaussian,
    BinnedClusterShearProfile,
    BinnedClusterNumberCounts,
    Likelihood,
    NamedParameters,
)



logger = logging.getLogger(__name__)


def _load_yaml_config(yml_file: str) -> Dict[str, Any]:
    """Load YAML file and return the CLPFirecrown section.

    Parameters
    ----------
    yml_file
        Path to YAML configuration file.

    Returns
    -------
    dict
        Parsed configuration dictionary (the CLPFirecrown sub-dictionary).
    """
    with open(yml_file, "r") as f:
        cfg = yaml.safe_load(f.read())
    # support configs that wrap settings under a top-level "CLPFirecrown" key
    return cfg.get("CLPFirecrown", cfg)


def _select_hmf(hmf_key: str, mass_def: str):
    """Return a pyccl halo mass function instance from a key."""
    hmf_dict = {
        "angulo12": ccl.halos.MassFuncAngulo12,
        "bocquet16": ccl.halos.MassFuncBocquet16,
        "bocquet20": ccl.halos.MassFuncBocquet20,
        "despali16": ccl.halos.MassFuncDespali16,
        "jenkins01": ccl.halos.MassFuncJenkins01,
        "press74": ccl.halos.MassFuncPress74,
        "sheth99": ccl.halos.MassFuncSheth99,
        "tinker08": ccl.halos.MassFuncTinker08,
        "tinker10": ccl.halos.MassFuncTinker10,
        "watson13": ccl.halos.MassFuncWatson13,
    }
    hmf_cls = hmf_dict.get(hmf_key, ccl.halos.MassFuncBocquet16)
    return hmf_cls(mass_def=mass_def)


def _build_completeness_and_purity(use_completeness, use_purity):
    """Create completeness and purity objects or None according to flags."""
    completeness = completeness_models.CompletenessAguena16() if use_completeness else None
    purity = purity_models.PurityAguena16() if use_purity else None
    return completeness, purity


def _populate_firecrown_parameters(yml_config: Dict[str, Any], set_params: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve Firecrown parameters: use sampled values from set_params when requested, otherwise fixed values."""
    firecrown_parameters = {}
    for param, value in yml_config.get("firecrown_parameters", {}).items():
        if value.get("sample", False):
            if param not in set_params:
                raise ValueError(f"Sampled Firecrown parameter '{param}' not provided in set_params")
            firecrown_parameters[param] = set_params[param]
        else:
            firecrown_parameters[param] = float(value.get("values"))
    return firecrown_parameters


def _build_mass_distribution(use_grid: bool, pivot_log_mass: float, pivot_z: float):
    """Return the appropriate mass proxy object depending on whether grid evaluation is chosen."""
    if use_grid:
        return mass_proxy.MurataUnbinned(pivot_log_mass=pivot_log_mass, pivot_redshift=pivot_z)
    return mass_proxy.MurataBinned(pivot_log_mass=pivot_log_mass, pivot_redshift=pivot_z)


def _apply_mass_distribution_parameters(mass_distribution, set_params: Dict[str, Any], firecrown_parameters: Dict[str, Any]):
    """Apply runtime-set attributes and firecrown parameters to the mass distribution object."""
    # set any attributes passed in set_params that exist on the object
    for k, v in set_params.items():
        if hasattr(mass_distribution, k):
            setattr(mass_distribution, k, v)

    # Fill mass-distribution parameter dictionary with defaults from firecrown parameters
    mass_distribution.parameters["mu0"] = firecrown_parameters.get("mass_distribution_mu0", 0.0)
    mass_distribution.parameters["mu1"] = firecrown_parameters.get("mass_distribution_mu1", 0.0)
    mass_distribution.parameters["mu2"] = firecrown_parameters.get("mass_distribution_mu2", 0.0)
    mass_distribution.parameters["sigma0"] = firecrown_parameters.get("mass_distribution_sigma0", 0.2)
    mass_distribution.parameters["sigma1"] = firecrown_parameters.get("mass_distribution_sigma1", 0.0)
    mass_distribution.parameters["sigma2"] = firecrown_parameters.get("mass_distribution_sigma2", 0.0)


def _apply_completeness_purity_parameters(completeness, purity, firecrown_parameters):
    """Apply firecrown-tuned parameters to completeness and purity objects when present."""
    if completeness is not None:
        completeness.parameters["a_n"] = firecrown_parameters.get("completeness_a_n", completeness.parameters["a_n"])
        completeness.parameters["b_n"] = firecrown_parameters.get("completeness_b_n", completeness.parameters["b_n"])
        completeness.parameters["a_logm_piv"] = firecrown_parameters.get("completeness_a_logm_piv", completeness.parameters["a_logm_piv"])
        completeness.parameters["b_logm_piv"] = firecrown_parameters.get("completeness_b_logm_piv", completeness.parameters["b_logm_piv"])
    if purity is not None:
        purity.parameters["a_n"] = firecrown_parameters.get("purity_a_n", purity.parameters["a_n"])
        purity.parameters["b_n"] = firecrown_parameters.get("purity_b_n", purity.parameters["b_n"])
        purity.parameters["a_logm_piv"] = firecrown_parameters.get("purity_a_logm_piv", purity.parameters["a_logm_piv"])
        purity.parameters["b_logm_piv"] = firecrown_parameters.get("purity_b_logm_piv", purity.parameters["b_logm_piv"])


def _make_binned_recipe(recipe_cls, cluster_theory, redshift_distribution, mass_distribution, completeness, purity,
                       mass_interval, true_z_interval, redshift_grid_size, mass_grid_size, proxy_grid_size):
    """Construct a binned cluster recipe (GridBinnedClusterRecipe or ExactBinnedClusterRecipe)."""
    return recipe_cls(
        cluster_theory=cluster_theory,
        redshift_distribution=redshift_distribution,
        mass_distribution=mass_distribution,
        completeness=completeness,
        purity=purity,
        mass_interval=mass_interval,
        true_z_interval=true_z_interval,
        redshift_grid_size=redshift_grid_size,
        mass_grid_size=mass_grid_size,
        proxy_grid_size=proxy_grid_size,
    )


def _build_shear_profile(cosmo_ccl, hmf, cluster_concentration, is_deltasigma, use_beta_interp, beta_parameters, set_params):
    """Construct and configure a ClusterShearProfile object."""
    shear = ClusterShearProfile(
        cosmo=cosmo_ccl,
        halo_mass_function=hmf,
        cluster_concentration=cluster_concentration,
        is_delta_sigma=is_deltasigma,
        use_beta_s_interp=use_beta_interp,
    )
    if "cluster_theory_cluster_concentration" in set_params:
        shear.cluster_concentration = set_params["cluster_theory_cluster_concentration"]
    if not is_deltasigma:
        shear.set_beta_parameters(*beta_parameters)
    if use_beta_interp:
        # The actual min/max z values will be set later when recipe is built; keep API analogous to old code
        # The calling code will set the interpolation domain on the recipe if needed.
        pass
    return shear


def build_cluster_recipes_from_config(yml_file: str, sacc_file: str, set_params: Dict[str, Any], **kwargs
                                     ) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Build cluster counts and (optionally) shear recipes from a YAML configuration and return theory/data vectors.

    This function:
    - Loads configuration from `yml_file` (supports a top-level CLPFirecrown key).
    - Builds a pyccl halo mass function and pyccl Cosmology (via build_ccl_cosmology_from_config).
    - Constructs completeness/purity models if requested.
    - Configures a mass-proxy model and fills its parameters from the `firecrown_parameters` block
      or from runtime `set_params`.
    - Builds either GridBinnedClusterRecipe or ExactBinnedClusterRecipe for counts and (optionally) for shear.
    - Reads the SACC file and returns theory predictions and data vectors.

    Parameters
    ----------
    yml_file
        Path to the YAML configuration.
    sacc_file
        Path to the SACC file used for data vectors.
    set_params
        Dictionary with runtime-supplied sampled parameters (e.g. from Firecrown/Cosmosis).
    **kwargs
        Currently unused; accepted for forward compatibility.

    Returns
    -------
    tuple
        (cluster_counts_theory, cluster_counts_data, cluster_shear_profile_theory, cluster_shear_data)
        cluster_shear_profile_theory and cluster_shear_data can be None if shear is not requested.
    """
    # Load config and resolve parameters
    yml_config = _load_yaml_config(yml_file)
    logger.debug("Loaded YAML config: %s", yml_config)

    # Basic parameters (with defaults)
    hmf_key = yml_config.get("hmf", "bocquet16")
    mass_def = str(yml_config.get("mass_def", "200c"))
    min_mass = yml_config.get("min_mass", 13.0)
    max_mass = yml_config.get("max_mass", 16.0)
    min_z = yml_config.get("min_z", 0.2)
    max_z = yml_config.get("max_z", 0.8)
    pivot_log_mass = yml_config.get("pivot_mass", 14.3)
    pivot_z = yml_config.get("pivot_z", 0.6)
    survey_name = yml_config.get("survey_name", "numcosmo_simulated_redshift_richness")

    use_shear_profile = yml_config.get("use_shear_profile", False)
    use_completeness = yml_config.get("use_completeness", None)
    use_purity = yml_config.get("use_purity", None)
    use_grid = yml_config.get("use_grid", True)
    is_deltasigma = yml_config.get("is_deltasigma", False)
    use_beta_interp = yml_config.get("use_beta_interp", False)
    redshift_grid_size = yml_config.get("redshift_grid_size", 20)
    mass_grid_size = yml_config.get("mass_grid_size", 60)
    proxy_grid_size = yml_config.get("proxy_grid_size", 20)
    beta_parameters = yml_config.get("beta_parameters", (10.0, 5.0))

    # Firecrown parameters (sampled or fixed)
    firecrown_parameters = _populate_firecrown_parameters(yml_config, set_params)
    cluster_concentration = firecrown_parameters.get("cluster_theory_cluster_concentration", None)

    # Build HMF and Cosmology
    hmf = _select_hmf(hmf_key, mass_def)
    cosmo_ccl = build_ccl_cosmology_from_config(yml_config, set_params)

    # completeness and purity
    completeness, purity = _build_completeness_and_purity(use_completeness, use_purity)

    # redshift distribution and abundance (cluster theory)
    redshift_distribution = kernel.SpectroscopicRedshift()
    abundance = ClusterAbundance(halo_mass_function=hmf, cosmo=cosmo_ccl)

    # mass distribution / proxy
    mass_distribution = _build_mass_distribution(use_grid, pivot_log_mass, pivot_z)
    _apply_mass_distribution_parameters(mass_distribution, set_params, firecrown_parameters)

    # apply completeness/purity firecrown parameters if present
    _apply_completeness_purity_parameters(completeness, purity, firecrown_parameters)

    # choose recipe class based on grid flag
    recipe_cls = GridBinnedClusterRecipe if use_grid else ExactBinnedClusterRecipe

    # build counts recipe
    counts_recipe = _make_binned_recipe(
        recipe_cls=recipe_cls,
        cluster_theory=abundance,
        redshift_distribution=redshift_distribution,
        mass_distribution=mass_distribution,
        completeness=completeness,
        purity=purity,
        mass_interval=(min_mass, max_mass),
        true_z_interval=(min_z, max_z),
        redshift_grid_size=redshift_grid_size,
        mass_grid_size=mass_grid_size,
        proxy_grid_size=proxy_grid_size,
    )

    # Build statistics object for counts and read SACC file
    average_on = ClusterProperty.NONE
    if yml_config.get("use_cluster_counts", True):
        average_on |= ClusterProperty.COUNTS
    if yml_config.get("use_mean_log_mass", False):
        average_on |= ClusterProperty.MASS
    if yml_config.get("use_mean_deltasigma", False):
        average_on |= ClusterProperty.DELTASIGMA
    if yml_config.get("use_mean_reduced_shear", False):
        average_on |= ClusterProperty.SHEAR

    counts_statistics = BinnedClusterNumberCounts(average_on, survey_name, counts_recipe)

    sacc_obj = sacc.Sacc.load_fits(sacc_file)
    counts_statistics.read(sacc_obj)

    logger.debug("Mass distribution parameters: %s", counts_statistics.cluster_recipe.mass_distribution.parameters)
    cluster_counts_theory = counts_statistics.get_binned_cluster_counts()
    cluster_counts_data = counts_statistics.data_vector

    # Optionally build shear recipe and statistics
    cluster_shear_profile_theory = None
    cluster_shear_data = None
    if use_shear_profile:
        shear = _build_shear_profile(
            cosmo_ccl=cosmo_ccl,
            hmf=hmf,
            cluster_concentration=cluster_concentration,
            is_deltasigma=is_deltasigma,
            use_beta_interp=use_beta_interp,
            beta_parameters=beta_parameters,
            set_params=set_params,
        )

        # If beta interpolation requires z range, set after knowing config's min/max
        if use_beta_interp and not is_deltasigma:
            # ClusterShearProfile API (in previous code) had set_beta_s_interp(min_z, max_z)
            try:
                shear.set_beta_s_interp(min_z, max_z)
            except Exception:
                # fall back silently if the method is not available or has different signature
                logger.debug("set_beta_s_interp not available or failed; continuing without setting interpolation domain.")

        # build shear recipe (purity intentionally None for shear as before)
        shear_recipe = _make_binned_recipe(
            recipe_cls=recipe_cls,
            cluster_theory=shear,
            redshift_distribution=redshift_distribution,
            mass_distribution=mass_distribution,
            completeness=completeness,
            purity=None,
            mass_interval=(min_mass, max_mass),
            true_z_interval=(min_z, max_z),
            redshift_grid_size=redshift_grid_size,
            mass_grid_size=mass_grid_size,
            proxy_grid_size=proxy_grid_size,
        )

        shear_statistics = BinnedClusterShearProfile(average_on, survey_name, shear_recipe)
        shear_statistics.read(sacc_obj)
        cluster_shear_profile_theory = shear_statistics.get_binned_cluster_property(average_on)
        cluster_shear_data = shear_statistics.data_vector

    return cluster_counts_theory, cluster_counts_data, cluster_shear_profile_theory, cluster_shear_data


def build_ccl_cosmology_from_config(yml_config: Dict[str, Any], set_params: Optional[Dict[str, Any]] = None) -> ccl.Cosmology:
    """
    Build a pyccl.Cosmology using config values + runtime sampled parameters.

    The configuration is expected to contain a "cosmological_parameters" block where each parameter is either
    marked as sampled (sample: true) in which case the value must be present in set_params, or fixed with a
    numeric "values" entry.

    Parameters
    ----------
    yml_config
        Configuration dictionary (the CLPFirecrown section).
    set_params
        Runtime-supplied sampled parameters (may be None).

    Returns
    -------
    pyccl.Cosmology
    """
    if set_params is None:
        set_params = {}

    # Start from a complete default cosmology
    cosmo_params = DEFAULT_CCL_PARAMS.copy()

    # Extract config cosmological parameters
    cosmological_block = yml_config.get("cosmological_parameters", {})

    for param, spec in cosmological_block.items():
        if spec.get("sample", False):
            # sampled → must come from set_params
            if param not in set_params:
                raise ValueError(f"Sampled cosmological parameter '{param}' not provided in set_params")
            cosmo_params[param] = float(set_params[param])
        else:
            # fixed → take from config
            cosmo_params[param] = float(spec["values"])

    # Filter only parameters pyccl actually accepts
    valid_ccl_params = {k: v for k, v in cosmo_params.items() if k in DEFAULT_CCL_PARAMS}

    return ccl.Cosmology(**valid_ccl_params)


# Default CCL parameter dictionary (used as base for building cosmologies)
DEFAULT_CCL_PARAMS = dict(
    Omega_c=0.2052,
    Omega_b=0.0448,
    h=0.71,
    n_s=0.963,
    sigma8=0.8,
    Omega_k=0.0,
    Neff=3.044,
    m_nu=0.0,
    w0=-1.0,
    wa=0.0,
    T_CMB=2.7255,
)
