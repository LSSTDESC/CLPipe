#!/usr/bin/env python
"""Convergence diagnostics: autocorrelation time (emcee).
"""
import sys

sys.path.append("../../")
import numpy as np
import emcee
from clpipe.utils.cosmosis_mcmc_plots import fits_to_walker_chain

CHAIN_PATH = "../../data/capish_simulation/chains/capish_chains_cc_wl.fits"
BURN_FRACTION = 0.15  # applied here at load time -- the saved fits is the FULL chain


def autocorr_report(chain, param_names):
    n_steps = chain.shape[0]
    rows = []
    for i, name in enumerate(param_names):
        sub = chain[:, :, i]
        try:
            tau = emcee.autocorr.integrated_time(sub)[0]
            converged = True
        except emcee.autocorr.AutocorrError:
            tau = emcee.autocorr.integrated_time(sub, tol=0)[0]
            converged = False
        ratio = n_steps / tau if tau and not np.isnan(tau) else np.nan
        rows.append((name, tau, ratio, converged))
    return rows


def main():
    chain, _, param_names = fits_to_walker_chain(CHAIN_PATH, burn_fraction=BURN_FRACTION)
    print("chain shape (n_steps, n_walkers, n_params), post-burn:", chain.shape)

    rows = autocorr_report(chain, param_names)

    print(f"{'param':>10s} {'tau':>8s} {'n_steps/tau':>12s} {'converged':>10s}")
    for name, tau, ratio, converged in rows:
        print(f"{name:>10s} {tau:8.1f} {ratio:12.1f} {str(converged):>10s}")

    valid_ratios = [r[2] for r in rows if not np.isnan(r[2])]
    worst_ratio = min(valid_ratios) if valid_ratios else float("nan")
    all_converged = all(r[3] for r in rows)
    print()
    print("OK" if all_converged else "!! CHECK convergence -- see module docstring !!")


if __name__ == "__main__":
    main()
