#!/usr/bin/env python
"""Ceci stage: fit NFW halo masses to CLMM cluster-ensemble profiles.
"""

import pickle
import warnings

import numpy as np

from ceci import PipelineStage
from ceci.file_types import PickleFile

from .file_types import FiducialCosmology


class CLClusterMass(PipelineStage):
    """Fit halo masses to stacked (and optionally individual) cluster profiles."""

    name = "CLClusterMass"

    inputs = [
        ("cluster_profiles", PickleFile),      # from CLClusterEnsembleProfiles
        ("fiducial_cosmology", FiducialCosmology),
    ]

    outputs = [
        ("cluster_masses", PickleFile),
    ]

    config_options = {
        # --- halo model / mass definition ---
        "halo_profile_model": "nfw",     # nfw, einasto, hernquist
        "delta": 200,                    # overdensity contrast
        "mass_def": "critical",          # "critical" or "mean"
        "concentration": 4.0,            # fixed c; if <= 0 use Bhattacharya13 c(M, z)
        # --- observable: DeltaSigma vs reduced shear ---
        "is_deltasigma": "auto",         # "auto" -> infer from profile_type; else "true"/"false"
        "z_src": 1.0,                    # effective source redshift (reduced shear only)
        # --- profile columns (must match CLClusterEnsembleProfiles) ---
        "tan_component": "tangential_comp",
        # --- covariance ---
        "cov_key": "auto",               # "auto" -> jackknife > bootstrap > sample; or tan_jk/tan_bs/tan_sc
        "use_covariance": True,          # full covariance vs diagonal-only
        # --- fitter (from clmm.support.sampler) ---
        "fitter": "curve_fit",           # "curve_fit" (fitters; gives errors),
                                         # "minimize" or "basinhopping" (samplers; peak only)
        "absolute_sigma": True,          # curve_fit only: sigma used in an absolute sense
        # --- radial range used in the fit (profile units; Mpc from upstream) ---
        "r_min_fit": 0.0,
        "r_max_fit": np.inf,
        # --- optimiser bounds on log10(M / Msun) ---
        "logm_min": 12.0,
        "logm_max": 16.0,
        "logm_0": 14.0,                  # initial guess (curve_fit p0 / sampler start)
        # --- individual cluster fits ---
        "fit_individual": False,
    }

    def run(self):
        self._setup_cosmology()

        # resolve model configuration once
        self.halo_profile_model = self.config["halo_profile_model"]
        self.delta = self.config["delta"]
        self.massdef = self.config["mass_def"].lower()
        self.fixed_concentration = self.config["concentration"]
        self.z_src = self.config["z_src"]
        self.tan_component = self.config["tan_component"]

        binned = pickle.load(open(self.get_input("cluster_profiles"), "rb"))

        results = {}
        for bin_key, entry in binned.items():
            ensemble = entry.get("clmm_cluster_ensemble")
            profile_type = entry.get("profile_type")
            is_ds = self._resolve_is_deltasigma(profile_type)

            out = {
                "cluster_bin_edges": entry.get("cluster_bin_edges"),
                "n_cl": entry.get("n_cl"),
                "profile_type": profile_type,
                "is_deltasigma": is_ds,
                "z_eff": None,
                "ensemble": None,
                "individual": None,
            }

            if ensemble is None:
                warnings.warn(f"Bin {bin_key}: no ensemble (n_cl <= 1), skipping fit.")
                results[bin_key] = out
                continue

            self.is_deltasigma = is_ds
            if not is_ds:
                warnings.warn(
                    f"Bin {bin_key}: fitting reduced shear with a single effective "
                    f"source redshift z_src={self.z_src}. Make sure this is set "
                    "appropriately; DeltaSigma is the primary supported path."
                )

            # effective (mean) cluster redshift of the bin
            z_eff = float(np.mean(self._col(ensemble.data, "z")))
            out["z_eff"] = z_eff

            # --- stacked / ensemble fit ---
            out["ensemble"] = self._fit_stacked(ensemble, z_eff)

            # --- individual fits (optional) ---
            if self.config["fit_individual"]:
                out["individual"] = self._fit_individual(ensemble)

            results[bin_key] = out
            self._log_bin(bin_key, out)

        pickle.dump(results, open(self.get_output("cluster_masses"), "wb"))

    # ------------------------------------------------------------------ #
    # setup helpers
    # ------------------------------------------------------------------ #
    def _setup_cosmology(self):
        """Load the fiducial cosmology as both a CCL and a CLMM object.

        The ceci-native ``FiducialCosmology`` wrapper (from ``.file_types``,
        shared with FirecrownPipeline) exposes the CCL-native parameters via
        ``.content``; older TXPipe-style wrappers instead provide ``.to_ccl()``.
        Both are handled here.
        """
        import clmm.cosmology.ccl

        with self.open_input("fiducial_cosmology", wrapper=True) as f:
            if hasattr(f, "to_ccl"):
                self.ccl_cosmo = f.to_ccl()
            else:
                self.ccl_cosmo = self._ccl_from_content(dict(f.content))

        self.clmm_cosmo = clmm.cosmology.ccl.CCLCosmology()
        self.clmm_cosmo.set_be_cosmo(self.ccl_cosmo)

    @staticmethod
    def _ccl_from_content(raw):
        """Build a ``pyccl.Cosmology`` from the fiducial-cosmology dict.

        Keys follow CCL-native naming (Omega_c, Omega_b, h, n_s, sigma8, and
        optionally Omega_k, w0, wa), the same set FirecrownPipeline reads.
        """
        import pyccl

        return pyccl.Cosmology(
            Omega_c=raw["Omega_c"],
            Omega_b=raw["Omega_b"],
            h=raw["h"],
            n_s=raw["n_s"],
            sigma8=raw["sigma8"],
            Omega_k=raw.get("Omega_k", 0.0),
            w0=raw.get("w0", -1.0),
            wa=raw.get("wa", 0.0),
        )

    def _resolve_is_deltasigma(self, profile_type):
        """Decide DeltaSigma (True) vs reduced shear (False).

        ``is_deltasigma = "auto"`` infers from the upstream ``profile_type``
        string; otherwise the config value ("true"/"false") wins.
        """
        cfg = self.config["is_deltasigma"]
        if isinstance(cfg, bool):
            return cfg
        if isinstance(cfg, str) and cfg.lower() != "auto":
            return cfg.lower() in ("true", "1", "yes", "deltasigma")
        pt = (profile_type or "").lower()
        return ("sigma" in pt) or ("excess" in pt)

    # ------------------------------------------------------------------ #
    # stacked (ensemble) fit
    # ------------------------------------------------------------------ #
    def _fit_stacked(self, ensemble, z_eff):
        radius = np.asarray(self._col(ensemble.stacked_data, "radius"), dtype=float)
        signal = np.asarray(
            self._col(ensemble.stacked_data, self.tan_component), dtype=float
        )
        cov, cov_key = self._get_covariance(ensemble, radius.size)

        # drop bins outside the fit range, with non-finite signal, or with a
        # non-positive covariance diagonal (e.g. empty bins kept upstream)
        mask = self._radial_mask(radius) & np.isfinite(signal) & (np.diag(cov) > 0)
        if mask.sum() < 2:
            warnings.warn("Fewer than 2 usable radial bins in stacked profile.")
            return None

        cov_f = cov[np.ix_(mask, mask)]
        if self.config["use_covariance"]:
            sigma = cov_f
        else:
            sigma = np.sqrt(np.clip(np.diag(cov_f), 1e-300, None))
        result = self._fit_mass(radius[mask], signal[mask], sigma, z_eff)
        result["cov_key"] = cov_key
        return result

    def _get_covariance(self, ensemble, nbins):
        """Select a computed covariance matrix from the ensemble.
        """
        cov_dict = getattr(ensemble, "cov", {}) or {}
        key = self.config["cov_key"]

        if key == "auto":
            key = None
            for candidate in ("tan_jk", "tan_bs", "tan_sc"):
                val = cov_dict.get(candidate)
                if val is not None and np.ndim(val) == 2:
                    key = candidate
                    break
        else:
            val = cov_dict.get(key)
            if val is None or np.ndim(val) != 2:
                key = None

        if key is None:
            shapes = {k: (None if v is None else np.shape(v)) for k, v in cov_dict.items()}
            warnings.warn(
                "No computed covariance on the ensemble (entries are None or not "
                "2-D); using identity covariance. The best-fit mass is then an "
                "unweighted fit and the reported error is NOT reliable. Set "
                "cov_type in CLClusterEnsembleProfiles (jackknife/bootstrap/sample) "
                f"to fix this. cov entries: {shapes}"
            )
            return np.eye(nbins), None

        return np.asarray(cov_dict[key], dtype=float), key

    # ------------------------------------------------------------------ #
    # individual cluster fits
    # ------------------------------------------------------------------ #
    def _fit_individual(self, ensemble):
        """Fit a mass to each individual cluster profile in the ensemble.

        Note
        ----
        IN CONSTRUCTION, NOT RELIABLE
        """
        data = ensemble.data
        n_cl = len(data)
        stacked_radius = np.asarray(self._col(ensemble.stacked_data, "radius"), dtype=float)
        stacked_cov, _ = self._get_covariance(ensemble, stacked_radius.size)
        approx_err = np.sqrt(np.clip(np.diag(stacked_cov), 0.0, None)) * np.sqrt(max(n_cl, 1))

        ids = self._col(data, "id", default=np.arange(n_cl))
        z_all = self._col(data, "z")
        out = []
        for i in range(n_cl):
            radius = np.asarray(data["radius"][i], dtype=float)
            signal = np.asarray(data[self.tan_component][i], dtype=float)
            z_cl = float(z_all[i])

            err = approx_err
            if err.shape[0] != radius.shape[0]:
                err = np.full(radius.shape[0], np.nanmedian(approx_err))

            mask = (
                self._radial_mask(radius)
                & np.isfinite(signal)
                & np.isfinite(err)
                & (err > 0)
            )
            if mask.sum() < 2:
                out.append({"cluster_id": _as_native(ids[i]), "z": z_cl, "logm": np.nan})
                continue

            res = self._fit_mass(radius[mask], signal[mask], err[mask], z_cl)
            res["cluster_id"] = _as_native(ids[i])
            out.append(res)
        return out

    # ------------------------------------------------------------------ #
    # fitting engine
    # ------------------------------------------------------------------ #
    def _fit_mass(self, radius, signal, sigma, z_cl):
        """Fit log10(mass) using the fitter selected by ``config['fitter']``.

        All three routines come from ``clmm.support.sampler``. ``sigma`` is a
        1-D array of per-bin errors or a 2-D covariance matrix.

        * ``curve_fit`` (``fitters``) returns ``(popt, pcov)``; the error is
          ``sqrt(pcov[0, 0])``. A 2-D ``sigma`` is used as the full covariance.
        * ``minimize`` / ``basinhopping`` (``samplers``) minimise the chi-square
          and return only the peak, so ``logm_err`` is NaN.
        """
        radius = np.asarray(radius, dtype=float)
        signal = np.asarray(signal, dtype=float)
        full_cov = np.ndim(sigma) == 2

        def model(r, log_m):
            return self._predict(r, log_m, z_cl)

        lo, hi = self.config["logm_min"], self.config["logm_max"]
        p0 = float(np.clip(self.config["logm_0"], lo, hi))
        fitter = str(self.config["fitter"]).lower()

        if fitter == "curve_fit":
            logm, logm_err, used_cov = self._fit_curve_fit(
                model, radius, signal, sigma, lo, hi, p0, full_cov
            )
        elif fitter in ("minimize", "basinhopping"):
            logm, logm_err, used_cov = self._fit_sampler(
                fitter, model, radius, signal, sigma, lo, hi, p0, full_cov
            )
        else:
            raise ValueError(
                f"Unknown fitter '{fitter}'. Options: curve_fit, minimize, basinhopping."
            )

        # chi-square at the best fit, for diagnostics
        res = signal - model(radius, logm)
        if used_cov:
            chi2 = float(res @ np.linalg.inv(sigma) @ res)
        else:
            err = np.sqrt(np.diag(sigma)) if full_cov else np.asarray(sigma, dtype=float)
            chi2 = float(np.sum((res / err) ** 2))
        dof = max(len(radius) - 1, 1)

        return {
            "logm": logm,
            "logm_err": logm_err,
            "m": 10 ** logm,
            "m_err": 10 ** logm * logm_err * np.log(10.0),
            "chi2": chi2,
            "dof": dof,
            "reduced_chi2": chi2 / dof,
            "concentration": float(self._concentration(logm, z_cl)),
            "z": float(z_cl),
            "n_radial_bins": int(len(radius)),
            "fitter": fitter,
            "used_covariance": bool(used_cov),
        }

    def _fit_curve_fit(self, model, radius, signal, sigma, lo, hi, p0, full_cov):
        """clmm.support.sampler.fitters['curve_fit'] -> (logm, logm_err, used_cov)."""
        from clmm.support.sampler import fitters

        abs_sigma = self.config["absolute_sigma"]
        try:
            popt, pcov = fitters["curve_fit"](
                model, radius, signal, sigma,
                absolute_sigma=abs_sigma, bounds=(lo, hi), p0=[p0],
            )
            used_cov = full_cov
        except Exception as exc:  # singular covariance / failed fit
            warnings.warn(
                f"curve_fit with the supplied errors failed ({exc}); "
                "retrying with diagonal errors."
            )
            diag = (
                np.sqrt(np.clip(np.diag(sigma), 1e-300, None))
                if full_cov
                else np.asarray(sigma, dtype=float)
            )
            popt, pcov = fitters["curve_fit"](
                model, radius, signal, diag,
                absolute_sigma=abs_sigma, bounds=(lo, hi), p0=[p0],
            )
            used_cov = False
        return float(popt[0]), float(np.sqrt(pcov[0][0])), used_cov

    def _fit_sampler(self, name, model, radius, signal, sigma, lo, hi, p0, full_cov):
        """clmm.support.sampler.samplers['minimize'|'basinhopping'].

        These minimise the chi-square objective and return only the best-fit
        point, so no parameter error is available (logm_err = NaN).
        """
        from clmm.support.sampler import samplers

        # build the chi-square weight (inverse covariance or inverse variance)
        used_cov = False
        if full_cov:
            try:
                weight = np.linalg.inv(sigma)
                used_cov = True
            except np.linalg.LinAlgError:
                warnings.warn("Singular covariance; using diagonal for the sampler.")
                weight = 1.0 / np.clip(np.diag(sigma), 1e-300, None)
        else:
            weight = 1.0 / np.asarray(sigma, dtype=float) ** 2

        def objective(params):
            log_m = float(np.ravel(params)[0])
            r = signal - model(radius, log_m)
            return float(r @ weight @ r) if used_cov else float(np.sum(r ** 2 * weight))

        if name == "minimize":
            kwargs = dict(method="L-BFGS-B", bounds=[(lo, hi)])
        else:  # basinhopping
            kwargs = dict(minimizer_kwargs=dict(method="L-BFGS-B", bounds=[(lo, hi)]))

        x = samplers[name](objective, np.array([p0]), **kwargs)
        warnings.warn(
            f"Fitter '{name}' returns only the best-fit point; logm_err is NaN. "
            "Use fitter='curve_fit' for a parameter error."
        )
        return float(np.ravel(x)[0]), float("nan"), used_cov

    # ------------------------------------------------------------------ #
    # model prediction
    # ------------------------------------------------------------------ #
    def _predict(self, radius, log_m, z_cl):
        """NFW DeltaSigma or reduced tangential shear at the given mass/redshift."""
        import clmm

        conc = self._concentration(log_m, z_cl)
        common = dict(
            r_proj=radius,
            mdelta=10 ** log_m,
            cdelta=conc,
            z_cl=z_cl,
            cosmo=self.clmm_cosmo,
            delta_mdef=self.delta,
            massdef=self.massdef,
            halo_profile_model=self.halo_profile_model,
        )
        if self.is_deltasigma:
            return clmm.compute_excess_surface_density(**common)
        return clmm.compute_reduced_tangential_shear(
            z_src=self.z_src, z_src_info="discrete", **common
        )

    def _concentration(self, log_m, z):
        """Fixed concentration, or Bhattacharya13 c(M, z) if ``concentration <= 0``."""
        if self.fixed_concentration is not None and self.fixed_concentration > 0:
            return self.fixed_concentration

        import pyccl

        rho_type = "matter" if self.massdef == "mean" else "critical"
        mass_def = pyccl.halos.MassDef(self.delta, rho_type)
        conc_model = pyccl.halos.concentration.ConcentrationBhattacharya13(mass_def=mass_def)
        a = 1.0 / (1.0 + z)
        try:
            return float(conc_model(self.ccl_cosmo, 10.0 ** log_m, a))
        except TypeError:
            # older pyccl API
            return float(conc_model._concentration(self.ccl_cosmo, 10.0 ** log_m, a))

    # ------------------------------------------------------------------ #
    # small utilities
    # ------------------------------------------------------------------ #
    def _radial_mask(self, radius):
        return (radius >= self.config["r_min_fit"]) & (radius <= self.config["r_max_fit"])

    @staticmethod
    def _col(table, name, default=None):
        """Column access with a couple of common fallbacks."""
        if name in table.colnames:
            return table[name]
        aliases = {
            "z": ("z_cl", "redshift"),
            "id": ("cluster_id", "unique_id"),
        }
        for alt in aliases.get(name, ()):
            if alt in table.colnames:
                return table[alt]
        if default is not None:
            return default
        raise KeyError(f"Column '{name}' not found (have: {table.colnames})")

    def _log_bin(self, bin_key, out):
        ens = out.get("ensemble")
        if ens is None:
            return
        print(
            f"[CLClusterMass] bin={bin_key} n_cl={out['n_cl']} "
            f"z_eff={out['z_eff']:.3f} "
            f"M={ens['m']:.3e} +/- {ens['m_err']:.3e} Msun "
            f"(logM={ens['logm']:.3f}+/-{ens['logm_err']:.3f}, "
            f"chi2/dof={ens['reduced_chi2']:.2f}, c={ens['concentration']:.2f})"
        )


def _as_native(value):
    """Convert numpy scalars to native Python types for a clean pickle."""
    try:
        return value.item()
    except AttributeError:
        return value
