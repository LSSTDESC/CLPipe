#!/usr/bin/env python
"""Regenerate ./chains/capish_chains_cc_wl.fits (used by analysis_summary.ipynb)
from the raw cosmosis chain.

Saves the FULL chain (burn_fraction=0.0) -- no burn-in is removed here.
Every consumer applies its own trim at load time instead: fits_to_samples()
and fits_to_walker_chain() both default to burn_fraction=0.15. Do not compute
statistics (plots, best-fit, autocorrelation, ...) directly off this file
without going through one of those two loaders.
"""
import numpy as np
import sys

sys.path.append("../../")
from clpipe.utils.cosmosis_mcmc_plots import chains_to_fits

DIR = ["./run_in2p3_both_mean/outputs_both/output_rp/number_counts_samples.txt"]

FIDUCIALS = [
    (r"\Omega_c", 0.22),
    (r"\sigma_8", 0.80),
    (r"\mu_0", 3.34),
    (r"\mu_m", 2.3 / np.log(10)),
    (r"\mu_z", 0.0),
    (r"\sigma_0", 0.56),
    (r"\sigma_m", 0.0),
    (r"\sigma_z", 0.0),
]


def main():
    chains_to_fits(
        paths=DIR,
        params=FIDUCIALS,
        labels=["capish_chains_cc_wl"],
        output_dir="../../data/capish_simulation/chains/",
        burn_fraction=0.0,
    )


if __name__ == "__main__":
    main()
