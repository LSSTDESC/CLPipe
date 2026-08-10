import ast
import configparser
import subprocess

import pytest
import yaml

from clpipe.clp_firecrown import CLPFirecrown
from .conftest import run_ceci_stage


COSMOLOGICAL_PARAMETERS = {
    "omega_c": {"sample": False, "values": 0.22},
    "sigma_8": {"sample": False, "values": 0.8},
    "omega_k": {"sample": False, "values": 0.0},
    "omega_b": {"sample": False, "values": 0.0448},
    "tau": {"sample": False, "values": 0.08},
    "n_s": {"sample": False, "values": 0.963},
    "h0": {"sample": False, "values": 0.71},
    "w": {"sample": False, "values": -1.0},
    "wa": {"sample": False, "values": 0.0},
}

FIRECROWN_PARAMETERS = {
    "mass_distribution_mu0": {"sample": True, "values": [2.0, 3.3439, 10.0]},
    "mass_distribution_mu1": {"sample": True, "values": [0.5, 0.9582, 2.0]},
    "mass_distribution_mu2": {"sample": True, "values": [-2.0, -0.0193, 2.0]},
    "mass_distribution_sigma0": {"sample": True, "values": [0.1, 0.5623, 2.0]},
    "mass_distribution_sigma1": {"sample": True, "values": [-0.6, 0.0455, 0.3]},
    "mass_distribution_sigma2": {"sample": True, "values": [-0.5, -0.0445, 2.0]},
    "completeness_a_n": {"sample": False, "values": 1.1321},
    "completeness_b_n": {"sample": False, "values": 0.7751},
    "completeness_a_logm_piv": {"sample": False, "values": 13.31},
    "completeness_b_logm_piv": {"sample": False, "values": 0.2025},
    "purity_a_n": {"sample": False, "values": 1.98298628209},
    "purity_b_n": {"sample": False, "values": 0.81212176229},
    "purity_a_logm_piv": {"sample": False, "values": 2.2183},
    "purity_b_logm_piv": {"sample": False, "values": -0.6592},
    "cluster_theory_cluster_concentration": {"sample": False, "values": 4},
}


def _base_firecrown_config(**overrides):
    cfg = {
        "hmf": "despali16",
        "min_mass": 12.0,
        "max_mass": 15.5,
        "min_z": 0.2,
        "max_z": 0.8,
        "mass_def": "200c",
        "use_shear_profile": True,
        "use_completeness": True,
        "use_purity": True,
        "use_grid": True,
        "is_deltasigma": True,
        "use_beta_interp": False,
        "beta_parameters": [10.0, 5.0],
        "pivot_mass": 14.3,
        "pivot_z": 0.5,
        "survey_name": "mock_survey",
        "sampler": "test",
        "filename": "output_rp/number_counts_samples.txt",
        "resume": False,
        "use_cluster_counts": True,
        "use_mean_log_mass": False,
        "use_mean_deltasigma": True,
        "emcee_walkers": 4,
        "emcee_samples": 8,
        "emcee_nsteps": 1,
        "cosmological_parameters": COSMOLOGICAL_PARAMETERS,
        "firecrown_parameters": FIRECROWN_PARAMETERS,
    }
    cfg.update(overrides)
    return cfg


def _make_stage(tmp_path, sacc_path, fiducial_cosmology_path, cfg):
    stage = CLPFirecrown(
        {
            "config": None,
            "clusters_sacc_file_cov": str(sacc_path),
            "fiducial_cosmology": str(fiducial_cosmology_path),
            "sampler_file": str(tmp_path / "sampler_file.ini"),
            "likelihood_file": str(tmp_path / "likelihood_file.py"),
            "priors_file": str(tmp_path / "priors_file.ini"),
        }
    )
    stage.config.update(cfg)
    return stage


def _without_param(params, name):
    d = dict(params)
    d.pop(name, None)
    return d


def _as_fixed(params):
    fixed = {}
    for name, spec in params.items():
        values = spec["values"]
        value = values[1] if isinstance(values, (list, tuple)) else values
        fixed[name] = {"sample": False, "values": value}
    return fixed


def _as_sampled(params, spread=0.1):
    sampled = {}
    for name, spec in params.items():
        values = spec["values"]
        if isinstance(values, (list, tuple)):
            sampled[name] = {"sample": True, "values": list(values)}
        else:
            delta = abs(values) * spread if values != 0 else spread
            sampled[name] = {"sample": True, "values": [values - delta, values, values + delta]}
    return sampled


def _section_sample_flags(ini_path, section, names):
    parser = configparser.ConfigParser()
    parser.read(ini_path)
    flags = {}
    for name in names:
        raw = parser[section][name]
        flags[name] = len(raw.split()) == 3
    return flags


def _write_stage_yaml(tmp_path, stage_name, cfg):
    config_path = tmp_path / "stage_config.yml"
    with open(config_path, "w") as f:
        yaml.safe_dump({stage_name: cfg}, f)
    return config_path


def _run_full_cosmosis(tmp_path, fiducial_cosmology_path, cfg, sacc_fixture):
    sacc_in_place = tmp_path / "clusters_sacc_file_cov.sacc"
    sacc_in_place.write_bytes(sacc_fixture.read_bytes())
    stage = _make_stage(tmp_path, sacc_in_place, fiducial_cosmology_path, cfg)
    stage.run()

    ini_path = tmp_path / "sampler_file.ini"
    return subprocess.run(
        ["cosmosis", str(ini_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=1800,
    )


# ---------------------------------------------------------------------------
# Tier 1: file generation only (ceci_stack)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "overrides",
    [
        {},
        {"use_shear_profile": False, "use_mean_deltasigma": False},
        {"use_cluster_counts": False},
        {"use_grid": False},
        {"is_deltasigma": False},
        {"use_completeness": False, "use_purity": False},
        {"use_beta_interp": True},
        {"two_halo_term": True, "boost_factor": True},
    ],
    ids=[
        "baseline", "counts_only", "shear_only", "exact_binned",
        "reduced_shear", "no_selection_function", "beta_interp",
        "two_halo_and_boost",
    ],
)
def test_generated_python_is_syntactically_valid(ceci_stack, tmp_path, mock_fiducial_cosmology, overrides):
    cfg = _base_firecrown_config(**overrides)
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    output_path = tmp_path / f"generated_{hash(frozenset(overrides.items())) & 0xffff}.py"
    ok = stage.generate_python_file(str(output_path))
    assert ok, f"generate_python_file failed for overrides={overrides}"

    source = output_path.read_text()
    tree = ast.parse(source)

    func_names = {
        node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
    }
    assert "build_likelihood" in func_names
    assert "get_cluster_recipe" in func_names

    if cfg["use_cluster_counts"]:
        assert "get_cluster_abundance" in func_names
        assert "BinnedClusterNumberCounts" in source
    if cfg["use_shear_profile"]:
        assert "get_cluster_shear_profile" in func_names
        assert "BinnedClusterShearProfile" in source
        assert "cluster_concentration=None" in source


def test_sacc_path_baked_in_is_basename_only(ceci_stack, tmp_path, mock_fiducial_cosmology):
    cfg = _base_firecrown_config()
    nested_sacc_path = tmp_path / "some" / "nested" / "dir" / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, nested_sacc_path, mock_fiducial_cosmology, cfg)

    output_path = tmp_path / "generated_basename_check.py"
    stage.generate_python_file(str(output_path))
    source = output_path.read_text()

    assert "sacc_path = 'clusters_sacc_file_cov.sacc'" in source
    assert str(nested_sacc_path) not in source
    assert "some" not in source and "nested" not in source


def test_full_matching_cosmological_parameters_still_generates(ceci_stack, tmp_path, mock_fiducial_cosmology):
    cfg = _base_firecrown_config()
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    params_path = tmp_path / "params_full_match.ini"
    assert stage.generate_cosmosis_parameters_file(str(params_path))


@pytest.mark.parametrize("name,bad_value", [("omega_c", 0.30), ("h0", 0.68), ("tau", 0.10)])
def test_cosmological_parameter_mismatch_cancels_generation(
    ceci_stack, tmp_path, mock_fiducial_cosmology, name, bad_value
):
    mismatched = dict(COSMOLOGICAL_PARAMETERS)
    mismatched[name] = {"sample": False, "values": bad_value}
    cfg = _base_firecrown_config(cosmological_parameters=mismatched)
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    params_path = tmp_path / f"params_mismatch_{name}.ini"
    assert not stage.generate_cosmosis_parameters_file(str(params_path))


# ---------------------------------------------------------------------------
# Tier 2: .ini content checks (full_stack, not slow)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("resume", [False, True])
def test_resume_flag_written_correctly(full_stack, tmp_path, mock_fiducial_cosmology, resume):
    cfg = _base_firecrown_config(resume=resume, filename=str(tmp_path / "chain.txt"))
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    ini_path = tmp_path / "test_resume.ini"
    assert stage.generate_ini_file(str(ini_path), "likelihood_file.py", "priors_file.ini")

    parser = configparser.ConfigParser()
    parser.read(ini_path)
    expected = "T" if resume else "F"
    assert parser["runtime"]["resume"] == expected


def test_filename_option_used_not_hardcoded_default(full_stack, tmp_path, mock_fiducial_cosmology):
    custom_filename = str(tmp_path / "my_custom_chain.txt")
    cfg = _base_firecrown_config(filename=custom_filename)
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    ini_path = tmp_path / "test_filename.ini"
    stage.generate_ini_file(str(ini_path), "likelihood_file.py", "priors_file.ini")

    parser = configparser.ConfigParser()
    parser.read(ini_path)
    assert parser["output"]["filename"] == custom_filename


def test_only_cosmo_sampled_writes_correct_priors(full_stack, tmp_path, mock_fiducial_cosmology):
    cfg = _base_firecrown_config(
        cosmological_parameters=_as_sampled(COSMOLOGICAL_PARAMETERS),
        firecrown_parameters=_as_fixed(FIRECROWN_PARAMETERS),
    )
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    params_path = tmp_path / "params_cosmo_only.ini"
    assert stage.generate_cosmosis_parameters_file(str(params_path))

    cosmo_flags = _section_sample_flags(
        params_path, "cosmological_parameters", COSMOLOGICAL_PARAMETERS.keys()
    )
    mor_flags = _section_sample_flags(
        params_path, "firecrown_number_counts", FIRECROWN_PARAMETERS.keys()
    )
    assert all(cosmo_flags.values())
    assert not any(mor_flags.values())


def test_only_mor_sampled_writes_correct_priors(full_stack, tmp_path, mock_fiducial_cosmology):
    cfg = _base_firecrown_config(
        cosmological_parameters=_as_fixed(COSMOLOGICAL_PARAMETERS),
        firecrown_parameters=_as_sampled(FIRECROWN_PARAMETERS),
    )
    sacc_path = tmp_path / "clusters_sacc_file_cov.sacc"
    stage = _make_stage(tmp_path, sacc_path, mock_fiducial_cosmology, cfg)

    params_path = tmp_path / "params_mor_only.ini"
    assert stage.generate_cosmosis_parameters_file(str(params_path))

    cosmo_flags = _section_sample_flags(
        params_path, "cosmological_parameters", COSMOLOGICAL_PARAMETERS.keys()
    )
    mor_flags = _section_sample_flags(
        params_path, "firecrown_number_counts", FIRECROWN_PARAMETERS.keys()
    )
    assert not any(cosmo_flags.values())
    assert all(mor_flags.values())


# ---------------------------------------------------------------------------
# Tier 3: stage.run() + actual cosmosis execution (full_stack, slow)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_generated_files_run_with_test_sampler(full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cfg = _base_firecrown_config(sampler="test", filename=str(tmp_path / "chain.txt"))
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"


@pytest.mark.slow
def test_only_cosmo_sampled_runs_successfully(full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        cosmological_parameters=_as_sampled(COSMOLOGICAL_PARAMETERS),
        firecrown_parameters=_as_fixed(FIRECROWN_PARAMETERS),
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"


@pytest.mark.slow
def test_only_mor_sampled_runs_successfully(full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense):
    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        cosmological_parameters=_as_fixed(COSMOLOGICAL_PARAMETERS),
        firecrown_parameters=_as_sampled(FIRECROWN_PARAMETERS),
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"


@pytest.mark.slow
@pytest.mark.parametrize(
    "missing_param",
    ["purity_a_n", "purity_b_n", "purity_a_logm_piv", "purity_b_logm_piv"],
)
def test_missing_purity_parameter_fails_at_runtime(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, missing_param
):
    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        use_purity=True,
        firecrown_parameters=_without_param(FIRECROWN_PARAMETERS, missing_param),
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    combined = result.stdout + result.stderr

    assert result.returncode != 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    assert missing_param in combined


@pytest.mark.slow
@pytest.mark.parametrize(
    "missing_param",
    [
        "completeness_a_n", "completeness_b_n",
        "completeness_a_logm_piv", "completeness_b_logm_piv",
    ],
)
def test_missing_completeness_parameter_fails_at_runtime(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, missing_param
):
    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        use_completeness=True,
        firecrown_parameters=_without_param(FIRECROWN_PARAMETERS, missing_param),
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    combined = result.stdout + result.stderr

    assert result.returncode != 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    assert missing_param in combined


@pytest.mark.slow
@pytest.mark.parametrize("is_deltasigma", [True, False], ids=["deltasigma", "reduced_shear"])
def test_missing_cluster_concentration_fails_at_runtime(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, is_deltasigma
):
    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        use_shear_profile=True,
        is_deltasigma=is_deltasigma,
        firecrown_parameters=_without_param(
            FIRECROWN_PARAMETERS, "cluster_theory_cluster_concentration"
        ),
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    combined = result.stdout + result.stderr

    assert result.returncode != 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    assert "cluster_theory_cluster_concentration" in combined


@pytest.mark.slow
def test_selection_function_disabled_runs_without_purity_completeness_params(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense
):
    stripped = FIRECROWN_PARAMETERS
    for name in [
        "purity_a_n", "purity_b_n", "purity_a_logm_piv", "purity_b_logm_piv",
        "completeness_a_n", "completeness_b_n",
        "completeness_a_logm_piv", "completeness_b_logm_piv",
    ]:
        stripped = _without_param(stripped, name)

    cfg = _base_firecrown_config(
        sampler="test",
        filename=str(tmp_path / "chain.txt"),
        use_purity=False,
        use_completeness=False,
        firecrown_parameters=stripped,
    )
    result = _run_full_cosmosis(tmp_path, mock_fiducial_cosmology, cfg, mock_cluster_sacc_dense)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"


# ---------------------------------------------------------------------------
# CLI dispatch
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_cli_invocation_generates_all_outputs(full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense):
    sacc_in_place = tmp_path / "clusters_sacc_file_cov.sacc"
    sacc_in_place.write_bytes(mock_cluster_sacc_dense.read_bytes())

    cfg = _base_firecrown_config(sampler="test", filename=str(tmp_path / "chain.txt"))
    config_path = _write_stage_yaml(tmp_path, "CLPFirecrown", cfg)

    ini_path = tmp_path / "sampler_file.ini"
    py_path = tmp_path / "likelihood_file.py"
    params_path = tmp_path / "priors_file.ini"

    result = run_ceci_stage(
        module="clpipe",
        stage_name="CLPFirecrown",
        config_path=config_path,
        io_args={
            "clusters_sacc_file_cov": str(sacc_in_place),
            "fiducial_cosmology": str(mock_fiducial_cosmology),
            "sampler_file": str(ini_path),
            "likelihood_file": str(py_path),
            "priors_file": str(params_path),
        },
        cwd=tmp_path,
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"

    assert ini_path.exists() and py_path.exists() and params_path.exists()
    ast.parse(py_path.read_text())

    parser = configparser.ConfigParser()
    parser.read(ini_path)
    assert "runtime" in parser