"""Tests for clpipe.cluster_theory_pred.

The main checks generate the likelihood and values files with CLPFirecrown
and verify that cluster_theory_pred builds the same model and uses the same
parameter values, so the two cannot drift apart.
"""
import configparser

import numpy as np
import pytest

# The module imports crow and firecrown at import time, so skip the whole
# file when the science stack is not installed (fast CI job).
pytest.importorskip("pyccl")
pytest.importorskip("crow")
pytest.importorskip("firecrown")

from firecrown.likelihood import NamedParameters, load_likelihood
from firecrown.updatable import ParamsMap

from clpipe.cluster_theory_pred import (
    build_cluster_recipes_from_config,
    compute_cluster_predictions,
)
from .test_clp_firecrown import (
    COSMOLOGICAL_PARAMETERS,
    FIRECROWN_PARAMETERS,
    _base_firecrown_config,
    _make_stage,
)


def _generated_theory(tmp_path, sacc_path, fiducial_cosmology_path, cfg, parameters, monkeypatch):
    """Theory vector and chi2 of the likelihood file written by CLPFirecrown."""
    stage = _make_stage(tmp_path, sacc_path, fiducial_cosmology_path, cfg)
    likelihood_path = tmp_path / "likelihood_file.py"
    assert stage.generate_python_file(str(likelihood_path))

    # Same build parameters as the [firecrown_likelihood] section of the
    # sampler file. The generated file opens the SACC file by basename.
    build_parameters = NamedParameters(
        {
            "use_cluster_counts": cfg["use_cluster_counts"],
            "use_mean_deltasigma": cfg["use_mean_deltasigma"],
            "use_mean_reduced_shear": not cfg["use_mean_deltasigma"],
            "use_mean_log_mass": cfg["use_mean_log_mass"],
        }
    )
    monkeypatch.chdir(sacc_path.parent)
    likelihood, tools = load_likelihood(str(likelihood_path), build_parameters)

    params_map = ParamsMap(parameters)
    likelihood.update(params_map)
    tools.update(params_map)
    tools.prepare()
    chi2 = likelihood.compute_chisq(tools)
    return likelihood.get_theory_vector(), chi2


@pytest.mark.slow
@pytest.mark.parametrize(
    "overrides",
    [
        {},
        {"use_shear_profile": False},
        {"use_purity": False, "use_completeness": False},
        {"hmf": "tinker08", "mass_grid_size": 40},
    ],
    ids=["counts_deltasigma", "counts_only", "no_selection", "tinker08"],
)
def test_prediction_matches_generated_likelihood(
    tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, monkeypatch, overrides
):
    cfg = _base_firecrown_config(**overrides)
    pred = compute_cluster_predictions(cfg, mock_cluster_sacc_dense, mock_fiducial_cosmology)

    expected, expected_chi2 = _generated_theory(
        tmp_path, mock_cluster_sacc_dense, mock_fiducial_cosmology, cfg, pred.parameters, monkeypatch
    )
    theory = [pred.counts.theory]
    if pred.shear is not None:
        theory.append(pred.shear.theory)
    np.testing.assert_allclose(np.concatenate(theory), expected, rtol=1e-12)
    assert pred.chi2 == pytest.approx(expected_chi2, rel=1e-12)
    assert (pred.shear is None) == (not cfg["use_shear_profile"])


def test_parameters_match_values_file(tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense):
    """The values CLPFirecrown writes for CosmoSIS are the ones used here."""
    cosmological_parameters = dict(COSMOLOGICAL_PARAMETERS)
    cosmological_parameters["sigma_8"] = {"sample": True, "values": [0.6, 0.83, 1.0]}
    cfg = _base_firecrown_config(cosmological_parameters=cosmological_parameters)

    pred = compute_cluster_predictions(cfg, mock_cluster_sacc_dense, mock_fiducial_cosmology)

    stage = _make_stage(tmp_path, mock_cluster_sacc_dense, mock_fiducial_cosmology, cfg)
    values_path = tmp_path / "priors_file.ini"
    assert stage.generate_cosmosis_parameters_file(str(values_path))
    values = configparser.ConfigParser()
    values.read(values_path)

    def start_value(raw):
        fields = raw.split()
        return float(fields[1] if len(fields) == 3 else fields[0])

    ccl_names = {"omega_c": "Omega_c", "omega_b": "Omega_b", "h0": "h", "n_s": "n_s",
                 "sigma_8": "sigma8", "omega_k": "Omega_k", "w": "w0", "wa": "wa"}
    for name, ccl_name in ccl_names.items():
        assert pred.parameters[ccl_name] == pytest.approx(
            start_value(values["cosmological_parameters"][name])
        ), name
    assert pred.parameters["sigma8"] == pytest.approx(0.83)
    for name in FIRECROWN_PARAMETERS:
        assert pred.parameters[name] == pytest.approx(
            start_value(values["firecrown_number_counts"][name])
        ), name


@pytest.mark.slow
def test_params_override_config_values(mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cfg = _base_firecrown_config(use_shear_profile=False)
    base = compute_cluster_predictions(cfg, mock_cluster_sacc_dense, mock_fiducial_cosmology)
    pred = compute_cluster_predictions(
        cfg,
        mock_cluster_sacc_dense,
        mock_fiducial_cosmology,
        params={
            "cosmological_parameters--sigma_8": 0.85,
            "Omega_c": 0.25,
            "firecrown_number_counts--mass_distribution_mu0": 3.5,
        },
    )
    assert pred.parameters["sigma8"] == 0.85
    assert pred.parameters["Omega_c"] == 0.25
    assert pred.parameters["mass_distribution_mu0"] == 3.5
    assert not np.allclose(pred.counts.theory, base.counts.theory)


def test_unknown_parameter_raises(mock_fiducial_cosmology, mock_cluster_sacc_dense):
    with pytest.raises(ValueError, match="Unknown parameter 'not_a_parameter'"):
        compute_cluster_predictions(
            _base_firecrown_config(),
            mock_cluster_sacc_dense,
            mock_fiducial_cosmology,
            params={"not_a_parameter": 1.0},
        )


def test_missing_cluster_parameter_raises(mock_fiducial_cosmology, mock_cluster_sacc_dense):
    firecrown_parameters = dict(FIRECROWN_PARAMETERS)
    firecrown_parameters.pop("purity_a_n")
    cfg = _base_firecrown_config(firecrown_parameters=firecrown_parameters)
    with pytest.raises(ValueError, match="purity_a_n"):
        compute_cluster_predictions(cfg, mock_cluster_sacc_dense, mock_fiducial_cosmology)


def test_fixed_cosmology_mismatch_raises(mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cosmological_parameters = dict(COSMOLOGICAL_PARAMETERS)
    cosmological_parameters["omega_c"] = {"sample": False, "values": 0.30}
    cfg = _base_firecrown_config(cosmological_parameters=cosmological_parameters)
    with pytest.raises(ValueError, match="does not match fiducial value"):
        compute_cluster_predictions(cfg, mock_cluster_sacc_dense, mock_fiducial_cosmology)


@pytest.mark.slow
def test_tuple_interface(mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cfg = _base_firecrown_config(use_shear_profile=False)
    counts_theory, counts_data, shear_theory, shear_data = build_cluster_recipes_from_config(
        cfg, mock_cluster_sacc_dense, fiducial_cosmology=mock_fiducial_cosmology
    )
    assert len(counts_theory) == len(counts_data) == 4
    assert shear_theory is None and shear_data is None
