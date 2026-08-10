"""
Shared fixtures for CLPipe tests.

All test artifacts (generated sacc/ini/py files) must live under pytest's
`tmp_path`, never under a path in the repo -- tmp_path is unique per test
and cleaned up automatically by pytest's own tmp dir rotation, so nothing
needs to be deleted manually.
"""
import subprocess
import sys

import numpy as np
import pytest
import sacc
import yaml
from scipy.linalg import block_diag


# ---------------------------------------------------------------------------
# Mock SACC builders
# ---------------------------------------------------------------------------

N_Z_BINS = 2
N_RICH_BINS = 2
N_RADIUS_BINS = 3
SURVEY_NAME = "mock_survey"
SURVEY_AREA = 20.0


COUNTS_PER_BIN = 40.0


def _radius_bin_names():
    return [f"radius_{i}" for i in range(N_RADIUS_BINS)]


def _correlated_block(n, variance, rho=0.3):
    """A simple constant-correlation covariance block: valid PSD for
    rho in (-1/(n-1), 1), which is satisfied for any rho in [0, 0.3] and
    n up to a few dozen. Used to fake a jackknife/bootstrap-like radius x
    radius correlation structure for a single (z, richness) stacked
    profile bin.
    """
    corr = np.full((n, n), rho)
    np.fill_diagonal(corr, 1.0)
    return variance * corr


def _build_mock_cluster_sacc(dense_shear_cov: bool, include_shear: bool = True) -> sacc.Sacc:
    """Build a mock cluster SACC file with cluster_counts (+ optionally
    cluster_delta_sigma) data points and a block-diagonal covariance.

    Args:
        dense_shear_cov: if True, the delta_sigma covariance block for each
            (z, richness) bin has genuine radius-radius correlation
            (constant-correlation model). If False, it's diagonal only
            (per-radius variance, no cross terms) -- mirrors the
            "diagonal_shear_covariance=True" CLPCovariance behavior.
        include_shear: if False, only cluster_counts data points are added
            (exercises the CLPCovariance `only_counts` branch).

    Returns:
        sacc.Sacc: populated, with covariance, NOT yet saved to disk.
    """
    s = sacc.Sacc()
    s.add_tracer("survey", SURVEY_NAME, SURVEY_AREA)

    z_edges = np.linspace(0.2, 0.8, N_Z_BINS + 1)
    rich_edges = np.linspace(1.0, 1.6, N_RICH_BINS + 1)  # log10(richness) edges
    radius_edges = np.geomspace(0.3, 3.0, N_RADIUS_BINS + 1)

    z_bins = [f"bin_z_{i}" for i in range(N_Z_BINS)]
    rich_bins = [f"bin_rich_{j}" for j in range(N_RICH_BINS)]
    radius_bins = _radius_bin_names()

    for i, name in enumerate(z_bins):
        s.add_tracer("bin_z", name, z_edges[i], z_edges[i + 1])
    for j, name in enumerate(rich_bins):
        s.add_tracer("bin_richness", name, rich_edges[j], rich_edges[j + 1])
    for k, name in enumerate(radius_bins):
        center = 0.5 * (radius_edges[k] + radius_edges[k + 1])
        s.add_tracer("bin_radius", name, radius_edges[k], radius_edges[k + 1], center)

    cluster_counts = sacc.standard_types.cluster_counts
    cluster_delta_sigma = sacc.standard_types.cluster_delta_sigma

    counts_values = []
    for z_name in z_bins:
        for r_name in rich_bins:
            s.add_data_point(cluster_counts, (SURVEY_NAME, r_name, z_name), COUNTS_PER_BIN)
            counts_values.append(COUNTS_PER_BIN)

    ds_blocks = []
    if include_shear:
        # A simple decaying delta-sigma profile, magnitude is not physically
        # calibrated -- these tests check covariance structure/propagation,
        # not lensing signal accuracy.
        for z_name in z_bins:
            for r_name in rich_bins:
                for k, rad_name in enumerate(radius_bins):
                    value = 50.0 / (k + 1)
                    s.add_data_point(
                        cluster_delta_sigma,
                        (SURVEY_NAME, r_name, z_name, rad_name),
                        value,
                    )
                variance = np.full(N_RADIUS_BINS, (0.05 * 50.0) ** 2)
                if dense_shear_cov:
                    ds_blocks.append(_correlated_block(N_RADIUS_BINS, variance[0], rho=0.3))
                else:
                    ds_blocks.append(np.diag(variance))

    counts_cov_block = np.diag(counts_values)
    if ds_blocks:
        full_cov = block_diag(counts_cov_block, *ds_blocks)
    else:
        full_cov = counts_cov_block

    s.add_covariance(full_cov)
    s.to_canonical_order()
    return s


@pytest.fixture
def mock_cluster_sacc_dense(tmp_path):
    """Mock cluster SACC (counts + delta_sigma) with a DENSE (correlated)
    delta_sigma covariance, written to tmp_path. Returns the file path.
    """
    s = _build_mock_cluster_sacc(dense_shear_cov=True, include_shear=True)
    path = tmp_path / "mock_clusters_dense.sacc"
    s.save_fits(str(path), overwrite=True)
    return path


@pytest.fixture
def mock_cluster_sacc_diagonal(tmp_path):
    """Mock cluster SACC (counts + delta_sigma) with a DIAGONAL-ONLY
    delta_sigma covariance, written to tmp_path. Returns the file path.
    """
    s = _build_mock_cluster_sacc(dense_shear_cov=False, include_shear=True)
    path = tmp_path / "mock_clusters_diagonal.sacc"
    s.save_fits(str(path), overwrite=True)
    return path


@pytest.fixture
def mock_cluster_sacc_counts_only(tmp_path):
    """Mock cluster SACC with ONLY cluster_counts (no shear/delta_sigma at
    all), written to tmp_path. Exercises the `only_counts` branch of
    CLPCovariance.run(), where merge_data_covariance should never be
    invoked (has_covariance and only_counts both matter -- here
    only_counts=True regardless of has_covariance).
    """
    s = _build_mock_cluster_sacc(dense_shear_cov=False, include_shear=False)
    path = tmp_path / "mock_clusters_counts_only.sacc"
    s.save_fits(str(path), overwrite=True)
    return path


# ---------------------------------------------------------------------------
# Minimal TJPCov / mor config fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mock_mor_parameters():
    return {
        "mass_func": "Despali16",
        "mass_def": "200c",
        "halo_bias": "Tinker10",
        "min_halo_mass": 1.0e12,
        "max_halo_mass": 3.16e15,
        "m_pivot": 14.3,
        "z_pivot": 0.5,
        "mu_p0": 3.3439,
        "mu_p1": 0.9582,
        "mu_p2": -0.0193,
        "sigma_p0": 0.5623,
        "sigma_p1": 0.0455,
        "sigma_p2": -0.0445,
    }


# ---------------------------------------------------------------------------
# Stack-availability fixtures
#
# Use these to gate tests by what they actually need, instead of skip logic
# scattered across classes:
#
# - ceci_stack: only `ceci` needs to be importable. Use for tests that just
#   construct a stage and generate output files (no tjpcov/crow/firecrown/
#   cosmosis execution). These should always run, even in a lightweight CI
#   job without the full DESC stack installed.
# - full_stack: the entire science stack needs to be importable. Use for
#   tests that actually execute covariance computation or run cosmosis.
#
# Pair full_stack tests with @pytest.mark.slow when they do real
# computation (covariance blocks, a cosmosis sampler run, etc.), so
# `pytest -m "not slow"` gives fast feedback even on a machine that has the
# full stack installed. Tests that only need full_stack to *import*
# something (e.g. generate_ini_file's `import cosmosis`/`import firecrown`)
# but don't execute it are not slow.
# ---------------------------------------------------------------------------

@pytest.fixture
def ceci_stack():
    pytest.importorskip("ceci")


@pytest.fixture
def full_stack():
    pytest.importorskip("ceci")
    pytest.importorskip("pyccl")
    pytest.importorskip("tjpcov")
    pytest.importorskip("crow")
    pytest.importorskip("firecrown")
    pytest.importorskip("cosmosis")


def run_ceci_stage(module, stage_name, config_path, io_args, cwd):
    """Run a single ceci PipelineStage via the CLI, matching how
    clpipe/__main__.py -> PipelineStage.main() is invoked in
    production.

    Args:
        module: python -m target, e.g. "clpipe"
        stage_name: e.g. "CLPCovariance"
        config_path: path to a YAML config file with a top-level key
            matching stage_name (same structure as your production configs)
        io_args: dict of {tag: path} for both inputs and outputs
        cwd: working directory to run the subprocess in (should be under
            tmp_path)

    Returns:
        subprocess.CompletedProcess
    """
    cmd = [sys.executable, "-m", module, stage_name, f"--config={config_path}"]
    for tag, path in io_args.items():
        cmd.append(f"--{tag}={path}")
    return subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=600
    )

@pytest.fixture
def mock_fiducial_cosmology(tmp_path):
    content = {
        "Omega_c": 0.22,
        "Omega_b": 0.0448,
        "h": 0.71,
        "n_s": 0.963,
        "sigma8": 0.8,
        "A_s": "nan",
        "Omega_k": 0.0,
        "Neff": 3.046,
        "w0": -1.0,
        "wa": 0.0,
    }
    path = tmp_path / "fiducial_cosmology.yml"
    with open(path, "w") as f:
        yaml.safe_dump(content, f)
    return path
