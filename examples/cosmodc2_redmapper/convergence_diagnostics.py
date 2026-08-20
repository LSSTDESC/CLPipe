#!/usr/bin/env python
"""Convergence diagnostics: autocorrelation time (emcee).

Loads the saved chain snapshots (../../data/cosmodc2_redmapper/chains/*.fits,
produced by save_chains.py) instead of the raw output_rp/number_counts_samples.txt
files.
"""
import sys

sys.path.append("../../")
import numpy as np
import emcee
from clpipe.utils.cosmosis_mcmc_plots import fits_to_walker_chain

CHAINS_DIR = "../../data/cosmodc2_redmapper/chains/"
BURN_FRACTION = 0.15  # applied here at load time -- the saved fits are the FULL chain

HMFS = ["bocquet16", "tinker08", "tinker10"]
RUNS = ["cosmo_0.22_0.7", "cosmo_0.22_0.9"]

mor_labels = ["mor_full", "mor_lensing", "mor_counts"] + [f"mor_{h}" for h in HMFS] + ["mor_exclude_first_richness_bin"]
both_labels = ["both_full", "both_lensing", "both_counts"] + [f"both_{h}" for h in HMFS] + [f"both_{run}" for run in RUNS] + ["both_exclude_first_richness_bin"]
cosmo_labels = ["cosmo_full", "cosmo_lensing", "cosmo_counts"]


def autocorr_report(chain, param_names):
    n_steps = chain.shape[0]
    rows = []
    for i, name in enumerate(param_names):
        sub = chain[:, :, i]
        try:
            tau = emcee.autocorr.integrated_time(sub)[0]  # tol=50 default -- raises if not yet reliable
            converged = True
        except emcee.autocorr.AutocorrError:
            tau = emcee.autocorr.integrated_time(sub, tol=0)[0]
            converged = False
        ratio = n_steps / tau if tau and not np.isnan(tau) else np.nan
        rows.append((name, tau, ratio, converged))
    return rows


def run_diagnostics(labels):
    for label in labels:
        path = f"{CHAINS_DIR}{label}.fits"
        try:
            chain, _, param_names = fits_to_walker_chain(path, burn_fraction=BURN_FRACTION)
        except Exception as e:
            print(f"{label:40s}  FAILED to load: {e}")
            continue

        rows = autocorr_report(chain, param_names)
        valid_ratios = [r[2] for r in rows if not np.isnan(r[2])]
        if not valid_ratios:
            print(f"{label:40s}  all tau estimates failed (chain too short/noisy)")
            continue
        worst_ratio = min(valid_ratios)
        all_converged = all(r[3] for r in rows)
        flag = "OK" if all_converged else "!! CHECK !!"
        print(f"{label:40s}  n_steps/tau(worst) = {worst_ratio:6.1f}   {flag}")

        if not all_converged:
            for name, tau, ratio, converged in rows:
                mark = "" if converged else "  <- not yet reliable per emcee (tol=50)"
                print(f"    {name:>10s}: tau={tau:8.1f}  n_steps/tau={ratio:6.1f}{mark}")


def main():
    print("=== BOTH (counts + lensing + MOR) ===")
    run_diagnostics(both_labels)

    print("\n=== MOR only ===")
    run_diagnostics(mor_labels)

    print("\n=== COSMO only ===")
    run_diagnostics(cosmo_labels)


if __name__ == "__main__":
    main()
