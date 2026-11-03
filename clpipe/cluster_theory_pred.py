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
from typing import Dict, Mapping, Optional, Tuple, Union

import numpy as np
from firecrown.connector.mapping import MappingCosmoSIS
from firecrown.likelihood import BinnedClusterShearProfile, NamedParameters, load_likelihood
from firecrown.updatable import ParamsMap


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
        matches = [(section, name) for section in values for name in values[section]
                   if key in (name, f"{section}--{name}")]
        if len(matches) != 1:
            raise ValueError(f"Parameter '{key}' does not match exactly one entry of the values file")
        section, name = matches[0]
        values[section][name] = float(value)


def _binned_prediction(likelihood, statistic) -> Dict[str, Optional[np.ndarray]]:
    """Prediction and measurement of one statistic, in the order of its SACC data vector.

    Keys: ``theory``, ``data``, ``covariance`` (block of the SACC
    covariance), ``errors`` (square root of its diagonal), ``z_edges`` and
    ``proxy_edges`` (log10 richness), both of shape (n, 2), and ``radius``
    (radius bin centers, None for counts).
    """
    bins = statistic.bins
    radius = None
    if isinstance(statistic, BinnedClusterShearProfile):
        radius = np.array([b.radius_center for b in bins])
    covariance = likelihood.get_cov(statistic)
    return {
        "theory": np.asarray(statistic.get_theory_vector(), dtype=float),
        "data": np.asarray(statistic.get_data_vector(), dtype=float),
        "covariance": covariance,
        "errors": np.sqrt(np.diag(covariance)),
        "z_edges": np.array([b.z_edges for b in bins]),
        "proxy_edges": np.array([b.mass_proxy_edges for b in bins]),
        "radius": radius,
    }


def compute_cluster_predictions(
    sampler_file: Union[str, os.PathLike],
    params: Optional[Mapping[str, float]] = None,
) -> Tuple[Optional[dict], Optional[dict], float, Dict[str, float]]:
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
    counts, shear
        Dictionaries described in :func:`_binned_prediction`, for the number
        counts (and mean log mass, if enabled) and the stacked shear profile.
        None when the likelihood does not have that statistic.
    chi2
        Chi-squared of the whole likelihood, using the full SACC covariance.
    parameters
        Values used, keyed ``section--name`` as in the CosmoSIS chain files.
    """
    run_dir = os.path.dirname(os.path.abspath(sampler_file))
    sampler = _read_ini(sampler_file)
    module = sampler["firecrown_likelihood"]  # options of the Firecrown CosmoSIS module

    values = _read_values(os.path.join(run_dir, sampler["pipeline"]["values"]))
    _apply_params(values, params or {})

    # From here on, the same steps as firecrown.connector.cosmosis.likelihood.
    # The generated build_likelihood only reads the use_* flags, written TRUE or FALSE.
    build_parameters = NamedParameters(
        {name: {"TRUE": True, "FALSE": False}.get(raw.upper(), raw) for name, raw in module.items()}
    )
    with contextlib.chdir(run_dir):  # the likelihood file opens the SACC file by basename
        likelihood, tools = load_likelihood(module["likelihood_source"], build_parameters)

    # Cluster parameters from the sampling sections, cosmology translated to
    # CCL names. omega_nu is set by the CosmoSIS consistency module in a chain.
    firecrown_params = {}
    for section in module["sampling_parameters_sections"].split():
        firecrown_params.update(values[section])
    mapping = MappingCosmoSIS()
    mapping.set_params_from_cosmosis(NamedParameters({"omega_nu": 0.0, **values["cosmological_parameters"]}))
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
    return counts, shear, float(chi2), parameters
