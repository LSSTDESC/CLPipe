"""Cluster theory predictions from the files written by CLPFirecrown.

``CLPFirecrown`` writes a Firecrown likelihood (``likelihood_file.py``), a
CosmoSIS values file (``priors_file.ini``) and a CosmoSIS ini file
(``sampler_file.ini``) that ties them together. This module loads that
likelihood, gives it parameter values the way the Firecrown CosmoSIS
connector does during a chain, and evaluates it at one point of parameter
space. Nothing about the model is rebuilt here, so the prediction follows
any change in ``CLPFirecrown`` or Firecrown.

Use it to compare the model with the data at the starting point of a chain,
at a chain best fit, or at any other set of parameters.

The only difference with a CosmoSIS run is the cosmology. Here CCL computes
the power spectrum and distances itself (with CAMB), while in the chain they
come from the CosmoSIS CAMB module, so the predictions can differ slightly
from the chain's.
"""

import configparser
import contextlib
import os
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Union

import numpy as np
from firecrown.connector.mapping import MappingCosmoSIS
from firecrown.likelihood import BinnedClusterShearProfile, NamedParameters, load_likelihood
from firecrown.modeling_tools import CCLCreationMode
from firecrown.updatable import ParamsMap

# Filled in by the CosmoSIS consistency module during a chain, but read by
# MappingCosmoSIS and not written to the values file by CLPFirecrown.
_COSMOSIS_DEFAULTS = {"omega_nu": 0.0}


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
        Number counts (and mean log mass, if enabled). None when the
        likelihood has no counts statistic.
    shear
        Stacked shear profile (DeltaSigma or reduced shear). None when the
        likelihood has no shear statistic.
    chi2
        Chi-squared of the whole likelihood, using the full SACC covariance.
    parameters
        Values used for the prediction, keyed ``section--name`` as in the
        CosmoSIS chain files.
    """

    counts: Optional[BinnedPrediction]
    shear: Optional[BinnedPrediction]
    chi2: float
    parameters: Dict[str, float]


def _ini_value(raw: str):
    """Convert a CosmoSIS ini entry to bool, int, float or str."""
    if raw.upper() in ("T", "TRUE", "Y", "YES"):
        return True
    if raw.upper() in ("F", "FALSE", "N", "NO"):
        return False
    for cast in (int, float):
        try:
            return cast(raw)
        except ValueError:
            pass
    return raw


def _read_ini(path: Union[str, os.PathLike]) -> configparser.ConfigParser:
    ini = configparser.ConfigParser(inline_comment_prefixes=(";", "#"), interpolation=None)
    if not ini.read(path):
        raise FileNotFoundError(path)
    return ini


def _read_values(path: Union[str, os.PathLike]) -> Dict[str, Dict[str, float]]:
    """Values file as {section: {name: value}}, with the starting value of sampled parameters."""
    ini = _read_ini(path)
    values = {}
    for section in ini.sections():
        values[section] = {}
        for name, raw in ini[section].items():
            fields = [float(field) for field in raw.split()]
            if len(fields) not in (1, 3):
                raise ValueError(f"Cannot read '{name} = {raw}' in [{section}] of {path}")
            values[section][name] = fields[0] if len(fields) == 1 else fields[1]
    return values


def _apply_params(values: Dict[str, Dict[str, float]], params: Mapping[str, float]) -> None:
    """Replace entries of ``values`` with ``params``, keyed ``section--name`` or ``name``."""
    for key, value in params.items():
        if "--" in key:
            section, name = key.split("--", 1)
            sections = [section] if name in values.get(section, {}) else []
        else:
            name = key
            sections = [section for section in values if name in values[section]]
        if len(sections) != 1:
            known = [f"{section}--{name}" for section in values for name in values[section]]
            reason = "is in several sections, add the section" if sections else "is not in the values file"
            raise ValueError(f"Parameter '{key}' {reason}. Known parameters: {', '.join(known)}")
        values[sections[0]][name] = float(value)


def _binned_prediction(likelihood, statistic) -> BinnedPrediction:
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
    sampler_file: Union[str, os.PathLike],
    params: Optional[Mapping[str, float]] = None,
) -> ClusterPredictions:
    """Evaluate the CLPFirecrown likelihood at one point and return predictions and data.

    Parameters
    ----------
    sampler_file
        CosmoSIS ini file written by CLPFirecrown (``sampler_file``). The
        likelihood file, values file and SACC file are read from the same
        folder, as when CosmoSIS runs there.
    params
        Values replacing the ones of the values file, for instance a chain
        best fit. Keys are CosmoSIS names, either ``section--name`` as in the
        chain files (``cosmological_parameters--sigma_8``) or just ``name``
        (``mass_distribution_mu0``). Without ``params``, the prediction is at
        the starting point of the chain.

    Returns
    -------
    ClusterPredictions
    """
    run_dir = os.path.dirname(os.path.abspath(sampler_file))
    sampler = _read_ini(sampler_file)
    if not sampler.has_section("firecrown_likelihood"):
        raise ValueError(f"{sampler_file} has no [firecrown_likelihood] section")
    module = sampler["firecrown_likelihood"]

    values = _read_values(os.path.join(run_dir, sampler["pipeline"]["values"]))
    _apply_params(values, params or {})

    # Same steps as firecrown.connector.cosmosis.likelihood: build the
    # likelihood with the module options, then collect the sampling sections
    # and, unless the likelihood makes its own cosmology, the CosmoSIS
    # cosmology translated to CCL names.
    build_parameters = NamedParameters({name: _ini_value(raw) for name, raw in module.items()})
    with contextlib.chdir(run_dir):  # the likelihood file opens the SACC file by basename
        likelihood, tools = load_likelihood(module["likelihood_source"], build_parameters)

    firecrown_params = {}
    for section in module.get("sampling_parameters_sections", "").split():
        firecrown_params.update(values[section])
    if tools.ccl_factory.creation_mode == CCLCreationMode.DEFAULT:
        mapping = MappingCosmoSIS()
        mapping.set_params_from_cosmosis(
            NamedParameters({**_COSMOSIS_DEFAULTS, **values["cosmological_parameters"]})
        )
        firecrown_params.update(mapping.asdict())
    params_map = ParamsMap(firecrown_params)
    params_map.use_lower_case_keys(True)

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
    parameters = {f"{section}--{name}": value for section in values for name, value in values[section].items()}
    return ClusterPredictions(counts=counts, shear=shear, chi2=float(chi2), parameters=parameters)
