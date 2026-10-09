"""Cluster theory predictions from a CLPFirecrown configuration.

This module builds the same Firecrown likelihood that ``CLPFirecrown``
writes to its ``likelihood_file`` output (same halo mass function,
mass-richness relation, selection function and CROW recipe), evaluates it
at one point of parameter space and returns the predicted and measured
data vectors. Use it to compare a configuration with the data at the
fiducial point, at a chain best fit, or at any other set of parameters.

The model options and their defaults come from ``CLPFirecrown`` itself, and
the parameter values follow the same rules as the CosmoSIS values file the
stage writes: the cosmology starts from the fiducial cosmology file, the
``cosmological_parameters`` block (CosmoSIS names) overrides it, and every
``firecrown_parameters`` entry takes its fixed value, or its starting value
when sampled.

The only difference with a CosmoSIS run is the linear power spectrum. Here
CCL computes it (``boltzmann_camb``), while in the chain it comes from the
CosmoSIS CAMB module, so the predictions can differ slightly from the
chain's.
"""

import math
import os
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Tuple, Union

import numpy as np
import pyccl as ccl
import sacc
import yaml
from crow import ClusterAbundance, ClusterShearProfile, kernel, mass_proxy
from crow import completeness_models, purity_models
from crow.properties import ClusterProperty
from crow.recipes.binned_exact import ExactBinnedClusterRecipe
from crow.recipes.binned_grid import GridBinnedClusterRecipe
from firecrown.likelihood import (
    BinnedClusterNumberCounts,
    BinnedClusterShearProfile,
    ConstGaussian,
)
from firecrown.modeling_tools import ModelingTools
from firecrown.updatable import ParamsMap

from .clp_firecrown import CLPFirecrown, _CCL_TO_COSMOSIS_COSMO_MAP

_COSMOSIS_TO_CCL_COSMO_MAP = {v: k for k, v in _CCL_TO_COSMOSIS_COSMO_MAP.items()}

_HALO_MASS_FUNCTIONS = {
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


@dataclass
class BinnedPrediction:
    """Prediction and measurement of one cluster statistic.

    All arrays follow the order of the statistic's data vector in the SACC
    file, so ``theory[i]``, ``data[i]`` and ``z_edges[i]`` refer to the same
    bin.

    Attributes
    ----------
    theory
        Predicted data vector.
    data
        Measured data vector.
    covariance
        Block of the SACC covariance for this statistic.
    z_edges
        Redshift bin edges, shape (n, 2).
    proxy_edges
        log10(richness) bin edges, shape (n, 2).
    radius
        Radius bin centers, shape (n,). None for the counts statistic.
    sacc_indices
        Indices of the data points in the SACC file.
    """

    theory: np.ndarray
    data: np.ndarray
    covariance: np.ndarray
    z_edges: np.ndarray
    proxy_edges: np.ndarray
    radius: Optional[np.ndarray]
    sacc_indices: np.ndarray

    @property
    def errors(self) -> np.ndarray:
        """Square root of the covariance diagonal."""
        return np.sqrt(np.diag(self.covariance))


@dataclass
class ClusterPredictions:
    """Output of :func:`compute_cluster_predictions`.

    Attributes
    ----------
    counts
        Number counts (and mean log mass, if enabled). None when
        ``use_cluster_counts`` is false.
    shear
        Stacked shear profile (DeltaSigma or reduced shear). None when
        ``use_shear_profile`` is false.
    chi2
        Chi-squared of the whole likelihood, using the full SACC covariance.
    parameters
        Parameter values used for the prediction, with CCL names for the
        cosmology and Firecrown names for the cluster parameters.
    """

    counts: Optional[BinnedPrediction]
    shear: Optional[BinnedPrediction]
    chi2: float
    parameters: Dict[str, float]


def _load_config(config: Union[str, os.PathLike, Mapping[str, Any]]) -> Dict[str, Any]:
    """Return the CLPFirecrown options, filling missing ones with the stage defaults.

    Parameters
    ----------
    config
        Path to a pipeline configuration file, or a dictionary. Both may hold
        the options under a top-level ``CLPFirecrown`` key.
    """
    if isinstance(config, Mapping):
        raw = dict(config)
    else:
        with open(config, "r") as f:
            raw = yaml.safe_load(f)
    section = raw.get("CLPFirecrown", raw)

    cfg = {name: option.default for name, option in CLPFirecrown.config_options.items()}
    cfg.update(section)
    return cfg


def _config_value(spec: Mapping[str, Any]) -> float:
    """Value of a parameter block entry: the starting value when sampled, else the fixed value."""
    values = spec["values"]
    if spec.get("sample", False):
        return float(values[1])
    return float(values)


def _halo_mass_function(cfg: Mapping[str, Any]):
    """Return the pyccl halo mass function selected by ``hmf``."""
    hmf_key = cfg["hmf"].lower()
    try:
        hmf_cls = _HALO_MASS_FUNCTIONS[hmf_key]
    except KeyError:
        raise ValueError(
            f"Unknown halo mass function '{hmf_key}'. "
            f"Available options are: {', '.join(sorted(_HALO_MASS_FUNCTIONS))}"
        ) from None
    return hmf_cls(mass_def=cfg["mass_def"])


def _cluster_recipe(cfg: Mapping[str, Any], cluster_theory, is_reduced_shear: bool = False):
    """Build the CROW recipe, as ``get_cluster_recipe`` in the generated likelihood file."""
    completeness = completeness_models.CompletenessAguena16() if cfg["use_completeness"] else None
    purity = purity_models.PurityAguena16LnProxy() if cfg["use_purity"] else None

    if is_reduced_shear:
        cluster_theory.set_beta_parameters(cfg["beta_parameters"][0], cfg["beta_parameters"][1])
        if cfg["use_beta_interp"]:
            cluster_theory.set_beta_s_interp(cfg["min_z"], cfg["max_z"])

    common = dict(
        cluster_theory=cluster_theory,
        redshift_distribution=kernel.SpectroscopicRedshift(),
        completeness=completeness,
        purity=purity,
        mass_interval=(cfg["min_mass"], cfg["max_mass"]),
        true_z_interval=(cfg["min_z"], cfg["max_z"]),
    )
    if cfg["use_grid"]:
        return GridBinnedClusterRecipe(
            mass_distribution=mass_proxy.MurataUnbinned(
                pivot_log_mass=cfg["pivot_mass"], pivot_redshift=cfg["pivot_z"]
            ),
            redshift_grid_size=cfg["redshift_grid_size"],
            mass_grid_size=cfg["mass_grid_size"],
            proxy_grid_size=cfg["proxy_grid_size"],
            **common,
        )
    return ExactBinnedClusterRecipe(
        mass_distribution=mass_proxy.MurataBinned(
            pivot_log_mass=cfg["pivot_mass"], pivot_redshift=cfg["pivot_z"]
        ),
        **common,
    )


def build_cluster_likelihood(
    config: Union[str, os.PathLike, Mapping[str, Any]], sacc_file: Union[str, os.PathLike]
) -> ConstGaussian:
    """Build the Firecrown likelihood that CLPFirecrown writes to ``likelihood_file``.

    The likelihood has read the SACC file but has not been updated with
    parameter values yet.

    Parameters
    ----------
    config
        Pipeline configuration file (or dictionary) with a ``CLPFirecrown``
        section, as used to run the stage.
    sacc_file
        SACC file with the cluster data vector and covariance
        (``clusters_sacc_file_cov``).

    Returns
    -------
    firecrown.likelihood.ConstGaussian
    """
    cfg = _load_config(config)
    if not (cfg["use_cluster_counts"] or cfg["use_shear_profile"]):
        raise ValueError("Set use_cluster_counts or use_shear_profile to build a likelihood.")

    # Same flags CLPFirecrown passes to build_likelihood through the sampler
    # file, where reduced shear is used whenever DeltaSigma is not.
    average_on = ClusterProperty.NONE
    if cfg["use_cluster_counts"]:
        average_on |= ClusterProperty.COUNTS
    if cfg["use_mean_log_mass"]:
        average_on |= ClusterProperty.MASS
    if cfg["use_mean_deltasigma"]:
        average_on |= ClusterProperty.DELTASIGMA
    else:
        average_on |= ClusterProperty.SHEAR

    # The placeholder cosmology is replaced by the ModelingTools one at
    # every evaluation.
    statistics = []
    if cfg["use_cluster_counts"]:
        abundance = ClusterAbundance(
            halo_mass_function=_halo_mass_function(cfg),
            cosmo=ccl.CosmologyVanillaLCDM(),
        )
        statistics.append(
            BinnedClusterNumberCounts(average_on, cfg["survey_name"], _cluster_recipe(cfg, abundance))
        )
    if cfg["use_shear_profile"]:
        shear_profile = ClusterShearProfile(
            cosmo=ccl.CosmologyVanillaLCDM(),
            halo_mass_function=_halo_mass_function(cfg),
            cluster_concentration=None,
            is_delta_sigma=cfg["is_deltasigma"],
            use_beta_s_interp=cfg["use_beta_interp"],
            two_halo_term=cfg["two_halo_term"],
            boost_factor=cfg["boost_factor"],
        )
        recipe = _cluster_recipe(cfg, shear_profile, is_reduced_shear=not cfg["is_deltasigma"])
        statistics.append(BinnedClusterShearProfile(average_on, cfg["survey_name"], recipe))

    likelihood = ConstGaussian(statistics)
    likelihood.read(sacc.Sacc.load_fits(str(sacc_file)))
    return likelihood


def _resolve_parameters(
    cfg: Mapping[str, Any],
    fiducial_cosmology: Union[str, os.PathLike],
    params: Optional[Mapping[str, float]],
    cosmology_defaults: Mapping[str, float],
    cluster_names,
) -> Dict[str, float]:
    """Parameter values for the prediction, as a flat dict for a Firecrown ParamsMap.

    Cosmology: Firecrown defaults < fiducial cosmology file < config
    ``cosmological_parameters`` < ``params``. Cluster parameters: config
    ``firecrown_parameters`` < ``params``.
    """
    cosmology = dict(cosmology_defaults)
    with open(fiducial_cosmology, "r") as f:
        fiducial = yaml.safe_load(f)
    cosmology.update({name: float(fiducial[name]) for name in cosmology if name in fiducial})

    for name, spec in cfg["cosmological_parameters"].items():
        if name == "tau":
            # CAMB input only, it does not enter the cluster prediction
            continue
        if name not in _COSMOSIS_TO_CCL_COSMO_MAP:
            raise ValueError(
                f"cosmological_parameters['{name}'] has no CCL equivalent. "
                f"Supported names are: {', '.join(sorted(_COSMOSIS_TO_CCL_COSMO_MAP))}"
            )
        ccl_name = _COSMOSIS_TO_CCL_COSMO_MAP[name]
        value = _config_value(spec)
        # Same check as CLPFirecrown.generate_cosmosis_parameters_file
        if not spec.get("sample", False) and not math.isclose(value, cosmology[ccl_name], rel_tol=1e-6):
            raise ValueError(
                f"cosmological_parameters['{name}']={spec['values']} "
                f"does not match fiducial value {cosmology[ccl_name]}"
            )
        cosmology[ccl_name] = value

    cluster = {name: _config_value(spec) for name, spec in cfg["firecrown_parameters"].items()}

    # Accept CosmoSIS names, optionally with their section prefix as in the
    # chain outputs (e.g. cosmological_parameters--omega_c), CCL names and
    # Firecrown names.
    for name, value in (params or {}).items():
        short = name.split("--")[-1]
        if short in _COSMOSIS_TO_CCL_COSMO_MAP:
            cosmology[_COSMOSIS_TO_CCL_COSMO_MAP[short]] = float(value)
        elif short in cosmology:
            cosmology[short] = float(value)
        elif short in cluster_names:
            cluster[short] = float(value)
        elif short != "tau":
            raise ValueError(f"Unknown parameter '{name}'")

    missing = sorted(set(cluster_names) - set(cluster))
    if missing:
        raise ValueError(
            f"No value for {', '.join(missing)}. Add them to firecrown_parameters "
            "in the configuration or pass them in params."
        )
    return {
        **{name: cosmology[name] for name in sorted(cosmology)},
        **{name: cluster[name] for name in sorted(cluster_names)},
    }


def _binned_prediction(likelihood: ConstGaussian, statistic) -> BinnedPrediction:
    bins = statistic.bins
    radius = None
    if isinstance(statistic, BinnedClusterShearProfile):
        radius = np.array([b.radius_center for b in bins])
    return BinnedPrediction(
        theory=np.asarray(statistic.get_theory_vector(), dtype=float),
        data=np.asarray(statistic.get_data_vector(), dtype=float),
        covariance=likelihood.get_cov(statistic),
        z_edges=np.array([b.z_edges for b in bins]),
        proxy_edges=np.array([b.mass_proxy_edges for b in bins]),
        radius=radius,
        sacc_indices=np.asarray(statistic.sacc_indices),
    )


def compute_cluster_predictions(
    config: Union[str, os.PathLike, Mapping[str, Any]],
    sacc_file: Union[str, os.PathLike],
    fiducial_cosmology: Union[str, os.PathLike],
    params: Optional[Mapping[str, float]] = None,
) -> ClusterPredictions:
    """Evaluate the CLPFirecrown likelihood at one point and return predictions and data.

    Parameters
    ----------
    config
        Pipeline configuration file (or dictionary) with a ``CLPFirecrown``
        section, as used to run the stage.
    sacc_file
        SACC file with the cluster data vector and covariance
        (``clusters_sacc_file_cov``).
    fiducial_cosmology
        Fiducial cosmology YAML file (CCL names), the ``fiducial_cosmology``
        input of the stage.
    params
        Values replacing the configuration ones, for instance a chain best
        fit. Keys can be CosmoSIS names (``omega_c``, ``sigma_8``, also with
        the ``cosmological_parameters--`` prefix of the chain files), CCL
        names (``Omega_c``, ``sigma8``) or Firecrown names
        (``mass_distribution_mu0``, also with the
        ``firecrown_number_counts--`` prefix).

    Returns
    -------
    ClusterPredictions
    """
    cfg = _load_config(config)
    likelihood = build_cluster_likelihood(cfg, sacc_file)
    tools = ModelingTools()

    values = _resolve_parameters(
        cfg,
        fiducial_cosmology,
        params,
        tools.required_parameters().get_default_values(),
        list(likelihood.required_parameters().get_params_names()),
    )
    params_map = ParamsMap(values)
    likelihood.update(params_map)
    tools.update(params_map)
    tools.prepare()
    chi2 = likelihood.compute_chisq(tools)

    counts = shear = None
    for guarded in likelihood.statistics:
        statistic = guarded.statistic
        if isinstance(statistic, BinnedClusterShearProfile):
            shear = _binned_prediction(likelihood, statistic)
        else:
            counts = _binned_prediction(likelihood, statistic)
    return ClusterPredictions(counts=counts, shear=shear, chi2=float(chi2), parameters=values)


def build_cluster_recipes_from_config(
    yml_file: Union[str, os.PathLike],
    sacc_file: Union[str, os.PathLike],
    set_params: Optional[Mapping[str, float]] = None,
    *,
    fiducial_cosmology: Union[str, os.PathLike],
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Return theory and data vectors as a tuple (older interface).

    Same as :func:`compute_cluster_predictions` with ``params=set_params``.

    Returns
    -------
    tuple
        (cluster_counts_theory, cluster_counts_data, cluster_shear_profile_theory, cluster_shear_data).
        The counts or shear entries are None when that statistic is not used.
    """
    pred = compute_cluster_predictions(yml_file, sacc_file, fiducial_cosmology, params=set_params)
    counts_theory = counts_data = shear_theory = shear_data = None
    if pred.counts is not None:
        counts_theory, counts_data = pred.counts.theory, pred.counts.data
    if pred.shear is not None:
        shear_theory, shear_data = pred.shear.theory, pred.shear.data
    return counts_theory, counts_data, shear_theory, shear_data
