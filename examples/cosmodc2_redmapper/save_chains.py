#!/usr/bin/env python
"""Regenerate every ./chains/*.fits snapshot used by plot_samples.ipynb,
summary_plot_best_fits.ipynb, and cov_matrices.ipynb from the raw cosmosis
chain files.
"""
import sys

sys.path.append("../../")
from clpipe.utils.cosmosis_mcmc_plots import chains_to_fits

DIR = "./baseline/"


def chain(subdir, run):
    return f"{DIR}/{subdir}/run_in2p3_{run}/outputs_{run}/output_rp/number_counts_samples.txt"


SUBDIRS = {
    "full": "cosmodc2_redmapper_full_analysis",
    "lensing": "cosmodc2_redmapper_lensing",
    "counts": "cosmodc2_redmapper_counts",
}

DIR_HMF = "./hmf_analysis/"


def chain_hmf(subdir, run):
    return f"{DIR_HMF}/{subdir}/run_in2p3_{run}/outputs_{run}/output_rp/number_counts_samples.txt"


HMFS = ["bocquet16", "tinker08", "tinker10", "watson13"]

FID_DIR = "./fiducial_cosmology_analysis"
RUNS = ["cosmo_0.22_0.7", "cosmo_0.22_0.9"]

EXCLUDE_BOTH = "./high_completeness_analysis/exclude_first_rich_bin/run_in2p3_both/outputs_both/output_rp/number_counts_samples.txt"
EXCLUDE_MOR = "./high_completeness_analysis/exclude_first_rich_bin/run_in2p3_mor/outputs_both/output_rp/number_counts_samples.txt"

PARAMS_MOR = [
    (r"\mu_0", 3.34),
    (r"\mu_m", 0.958236982),
    (r"\mu_z", -0.0192802),
    (r"\sigma_0", 0.562317194),
    (r"\sigma_m", 0.04552506),
    (r"\sigma_z", -0.0445),
]
PARAMS_BOTH = [(r"\Omega_c", 0.22), (r"\sigma_8", 0.80), *PARAMS_MOR]
PARAMS_COSMO = [(r"\Omega_c", 0.22), (r"\sigma_8", 0.80)]


def main():
    mor_paths = (
        [chain(SUBDIRS[k], "mor") for k in ("full", "lensing", "counts")]
        + [chain_hmf(h, "mor") for h in HMFS]
        + [EXCLUDE_MOR]
    )
    mor_labels = (
        ["mor_full", "mor_lensing", "mor_counts"]
        + [f"mor_{h}" for h in HMFS]
        + ["mor_exclude_first_richness_bin"]
    )

    both_paths = (
        [chain(SUBDIRS[k], "both") for k in ("full", "lensing", "counts")]
        + [chain_hmf(h, "both") for h in HMFS]
        + [f"{FID_DIR}/{run}/run_in2p3_both/outputs_both/output_rp/number_counts_samples.txt" for run in RUNS]
        + [EXCLUDE_BOTH]
    )
    both_labels = (
        ["both_full", "both_lensing", "both_counts"]
        + [f"both_{h}" for h in HMFS]
        + [f"both_{run}" for run in RUNS]
        + ["both_exclude_first_richness_bin"]
    )

    cosmo_paths = [chain(SUBDIRS[k], "cosmo") for k in ("full", "lensing", "counts")]
    cosmo_labels = ["cosmo_full", "cosmo_lensing", "cosmo_counts"]

    chains_to_fits(paths=mor_paths, params=PARAMS_MOR, labels=mor_labels, output_dir="../../data/cosmodc2_redmapper/chains/", burn_fraction=0.0)
    chains_to_fits(paths=both_paths, params=PARAMS_BOTH, labels=both_labels, output_dir="../../data/cosmodc2_redmapper/chains/", burn_fraction=0.0)
    chains_to_fits(paths=cosmo_paths, params=PARAMS_COSMO, labels=cosmo_labels, output_dir="../../data/cosmodc2_redmapper/chains/", burn_fraction=0.0)


if __name__ == "__main__":
    main()
