import numpy as np
import pytest
import sacc
import yaml

from clpipe.clp_covariance import CLPCovariance
from .conftest import run_ceci_stage

pytestmark = pytest.mark.filterwarnings(
    "ignore:Empty index selected.*:UserWarning"
)


def _base_tjpcov_config(mor_parameters, **overrides):
    cfg = {
        "replace_tjpcov_cov": True,
        "sel_func": True,
        "wazp_catalog": False,
        "diagonal_shear_covariance": True,
        "use_mpi": False,
        "do_xi": False,
        "cov_type": ["ClusterCountsGaussian", "ClusterCountsSSC"],
        "photo-z": {"sigma_0": 0.05},
        "mor_parameters": mor_parameters,
    }
    cfg.update(overrides)
    return cfg


def _run_tjpcov_stage(input_sacc_path, output_sacc_path, fiducial_cosmology_path, tjpcov_cfg):
    stage = CLPCovariance(
        {
            "config": None,
            "clusters_sacc_file": str(input_sacc_path),
            "clusters_sacc_file_cov": str(output_sacc_path),
            "fiducial_cosmology": str(fiducial_cosmology_path),
        }
    )
    stage.config.update(tjpcov_cfg)
    stage.run()
    return sacc.Sacc.load_fits(str(output_sacc_path))


def _write_stage_yaml(tmp_path, stage_name, cfg):
    config_path = tmp_path / "stage_config.yml"
    with open(config_path, "w") as f:
        yaml.safe_dump({stage_name: cfg}, f)
    return config_path


def _diag_block(cov, idx):
    idx = np.asarray(sorted(idx))
    return cov[np.ix_(idx, idx)]


def _radius_group_indices(sacc_obj, data_type):
    data_points = sacc_obj.get_data_points(data_type=data_type)
    groups = {}
    for dp in data_points:
        key = dp.tracers[:3]
        idx = sacc_obj.indices(data_type=data_type, tracers=dp.tracers)[0]
        groups.setdefault(key, set()).add(idx)
    return groups


@pytest.mark.slow
def test_dense_merge_preserves_off_diagonal(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    output_path = tmp_path / "cov_dense_kept.sacc"
    cfg = _base_tjpcov_config(mock_mor_parameters, diagonal_shear_covariance=False)
    out = _run_tjpcov_stage(mock_cluster_sacc_dense, output_path, mock_fiducial_cosmology, cfg)

    cds = sacc.standard_types.cluster_delta_sigma
    groups = _radius_group_indices(out, cds)
    idx = next(iter(groups.values()))
    block = _diag_block(out.covariance.covmat, idx)
    off_diag_max = np.max(np.abs(block - np.diag(np.diag(block))))
    assert off_diag_max > 0, "expected preserved radius-radius correlation"


@pytest.mark.slow
def test_diagonal_only_merge_drops_off_diagonal(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    output_path = tmp_path / "cov_diag_only.sacc"
    cfg = _base_tjpcov_config(mock_mor_parameters, diagonal_shear_covariance=True)
    out = _run_tjpcov_stage(mock_cluster_sacc_dense, output_path, mock_fiducial_cosmology, cfg)

    cds = sacc.standard_types.cluster_delta_sigma
    groups = _radius_group_indices(out, cds)
    idx = next(iter(groups.values()))
    block = _diag_block(out.covariance.covmat, idx)
    off_diag_max = np.max(np.abs(block - np.diag(np.diag(block))))
    assert off_diag_max == 0.0, "expected off-diagonal terms stripped"
    assert np.all(np.diag(block) > 0), "variances should still be kept"


@pytest.mark.slow
def test_counts_only_input_completes_without_merge(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_counts_only, mock_mor_parameters
):
    output_path = tmp_path / "cov_counts_only.sacc"
    cfg = _base_tjpcov_config(mock_mor_parameters)
    out = _run_tjpcov_stage(mock_cluster_sacc_counts_only, output_path, mock_fiducial_cosmology, cfg)

    assert out.get_data_types() == [sacc.standard_types.cluster_counts]
    assert out.covariance is not None


@pytest.mark.slow
def test_mixed_input_data_types_preserved(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    output_path = tmp_path / "cov_mixed.sacc"
    cfg = _base_tjpcov_config(mock_mor_parameters)
    out = _run_tjpcov_stage(mock_cluster_sacc_dense, output_path, mock_fiducial_cosmology, cfg)

    assert set(out.get_data_types()) == {
        sacc.standard_types.cluster_counts,
        sacc.standard_types.cluster_delta_sigma,
    }


@pytest.mark.slow
def test_ssc_replaced_counts_variance_at_least_poisson(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    output_path = tmp_path / "cov_ssc_check.sacc"
    cfg = _base_tjpcov_config(mock_mor_parameters, replace_tjpcov_cov=True)
    out = _run_tjpcov_stage(mock_cluster_sacc_dense, output_path, mock_fiducial_cosmology, cfg)

    cc = sacc.standard_types.cluster_counts
    counts_points = out.get_data_points(data_type=cc)

    for dp in counts_points:
        idx = out.indices(data_type=cc, tracers=dp.tracers)[0]
        replaced_variance = out.covariance.covmat[idx, idx]
        poisson_variance = dp.value
        assert replaced_variance >= poisson_variance, (
            f"SSC-replaced variance ({replaced_variance}) < Poisson "
            f"variance ({poisson_variance}) for tracers {dp.tracers}."
        )


@pytest.mark.slow
def test_diagonal_shear_covariance_is_noop_for_counts_only(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_counts_only, mock_mor_parameters
):
    cfg_true = _base_tjpcov_config(mock_mor_parameters, diagonal_shear_covariance=True)
    cfg_false = _base_tjpcov_config(mock_mor_parameters, diagonal_shear_covariance=False)
    out_true = _run_tjpcov_stage(
        mock_cluster_sacc_counts_only, tmp_path / "cov_true.sacc", mock_fiducial_cosmology, cfg_true
    )
    out_false = _run_tjpcov_stage(
        mock_cluster_sacc_counts_only, tmp_path / "cov_false.sacc", mock_fiducial_cosmology, cfg_false
    )
    np.testing.assert_allclose(
        out_true.covariance.covmat, out_false.covariance.covmat
    )


@pytest.mark.slow
def test_empty_index_warning_only_from_cluster_counts_lookup(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    cfg = _base_tjpcov_config(mock_mor_parameters)
    output_path = tmp_path / "cov_warning_check.sacc"
    out = _run_tjpcov_stage(mock_cluster_sacc_dense, output_path, mock_fiducial_cosmology, cfg)

    cds = sacc.standard_types.cluster_delta_sigma
    groups = _radius_group_indices(out, cds)
    assert all(len(idx) > 0 for idx in groups.values()), (
        "cluster_delta_sigma tracer group came back empty"
    )


@pytest.mark.slow
def test_missing_fiducial_cosmology_input_fails_clearly(full_stack, tmp_path, mock_cluster_sacc_dense, mock_mor_parameters):
    cfg = _base_tjpcov_config(mock_mor_parameters)
    output_path = tmp_path / "cov_no_fiducial.sacc"
    with pytest.raises(Exception):
        CLPCovariance(
            {
                "config": None,
                "clusters_sacc_file": str(mock_cluster_sacc_dense),
                "clusters_sacc_file_cov": str(output_path),
            }
        )


# ---------------------------------------------------------------------------
# CLI dispatch
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_cli_invocation_produces_valid_covariance(
    full_stack, tmp_path, mock_fiducial_cosmology, mock_cluster_sacc_dense, mock_mor_parameters
):
    cfg = _base_tjpcov_config(mock_mor_parameters)
    config_path = _write_stage_yaml(tmp_path, "CLPCovariance", cfg)
    output_path = tmp_path / "cov_cli.sacc"

    result = run_ceci_stage(
        module="clpipe",
        stage_name="CLPCovariance",
        config_path=config_path,
        io_args={
            "clusters_sacc_file": str(mock_cluster_sacc_dense),
            "clusters_sacc_file_cov": str(output_path),
            "fiducial_cosmology": str(mock_fiducial_cosmology),
        },
        cwd=tmp_path,
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"

    out = sacc.Sacc.load_fits(str(output_path))
    assert out.covariance is not None