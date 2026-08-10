"""
MCMC triangle plots from CosmoSIS postprocess chains,
plus a helper to persist chains as FITS tables for later use.

Typical usage
─────────────
from triangle_plot import ChainStyle, PlotConfig, plot_triangle, chains_to_fits

# ── Plot ──────────────────────────────────────────────────────────────────────
styles = [
    ChainStyle("Counts + Lensing", filled=True),
    ChainStyle("Lensing"),
    ChainStyle("Cluster Counts"),
]

samples = plot_triangle(
    paths   = [path_full, path_lensing, path_counts],
    params  = PARAMS_MOR,          # list of (name, fiducial_value)
    styles  = styles,
    name    = "mor_baseline",
    config  = PlotConfig(figsize=7, output_dir="./figures"),
)

# ── Save chains to FITS (from CosmoSIS paths) ─────────────────────────────────
chains_to_fits(
    paths      = [path_full, path_lensing, path_counts],
    params     = PARAMS_MOR,
    labels     = ["full", "lensing", "counts"],
    output_dir = "./chains",
)

# ── Or save MCSamples you already loaded (e.g. after plot_triangle) ───────────
chains_to_fits(
    samples    = samples,
    params     = PARAMS_MOR,
    labels     = ["full", "lensing", "counts"],
    output_dir = "./chains",
)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from getdist import MCSamples, plots
from IPython.display import Math, display

try:
    from cosmosis.postprocessing.inputs import read_input as _cosmosis_read
except ImportError:  # pragma: no cover
    _cosmosis_read = None

DEFAULT_COLORS = ['#000000', '#e41a1c', '#377eb8']

GETDIST_SETTINGS = {
    "mult_bias_correction_order": 1,
    "smooth_scale_2D": 4,
    "smooth_scale_1D": 3,
    "boundary_correction_order": 1,
}

# RC overrides applied only while the function runs (restored afterwards)
_PAPER_RC = {
    "font.family":       "serif",
    "font.serif":        ["Computer Modern Roman", "DejaVu Serif"],
    "text.usetex":       False,   # flip to True if LaTeX is available
    "axes.labelsize":    18,
    "xtick.labelsize":    13,
    "ytick.labelsize":    13,
    "legend.fontsize":   15,
    "legend.frameon":   False,
    "figure.dpi":       150,
}


@dataclass
class ChainStyle:
    """
    Visual style for one chain in the triangle plot.

    Parameters
    ----------
    label       : Legend label for this chain.
    filled      : Fill the 68 % / 95 % contours (good for the "combined" chain).
    fill_alpha  : Opacity of the filled region (ignored when filled=False).
    line_alpha  : Opacity of the contour lines.
    linestyle   : '-' solid  |  '--' dashed  |  ':' dotted.
    linewidth   : Contour line width.
    color       : Hex color.  None → taken from DEFAULT_COLORS in order.
    """
    label:      str
    filled:     bool  = False
    fill_alpha: float = 0.25
    line_alpha: float = 0.90
    linestyle:  str   = "-"
    linewidth:  float = 1.8
    color:      Optional[str] = None


@dataclass
class PlotConfig:
    """
    Global knobs for the triangle plot.

    Parameters
    ----------
    figsize         : Side length in inches (square).  7 ≈ 6-param, 8 ≈ 8-param.
    burn_fraction   : Fraction of each chain to discard as burn-in.
    param_limits    : GetDist-style dict, e.g. {r'\\Omega_c': [0.1, 0.4]}.
    output_dir      : Where to save the PDF.  Defaults to current directory.
    dpi             : DPI for on-screen / PNG preview (PDF is always vector).
    print_constraints : Print LaTeX mean ± σ for every parameter after plotting.
    """
    figsize:           float       = 7.0
    burn_fraction:     float       = 0.3
    param_limits:      dict        = field(default_factory=dict)
    output_dir:        str         = "."
    dpi:               int         = 150
    print_constraints: bool        = True


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

def _load_cosmosis_chain(
    path: str,
    param_names: list[str],
    burn_fraction: float,
) -> MCSamples:
    """Read one CosmoSIS chain file and return a trimmed MCSamples object."""
    if _cosmosis_read is None:
        raise ImportError("CosmoSIS is not installed; cannot read chain files.")

    raw   = _cosmosis_read(path)
    rows  = raw[1]["data"][1][0]
    n_par = len(param_names)
    array = np.array([[row[i] for i in range(n_par)] for row in rows])

    samples = MCSamples(
        samples  = array,
        names    = param_names,
        labels   = param_names,
        settings = GETDIST_SETTINGS,
    )
    samples.removeBurn(burn_fraction)
    return samples


def _resolve_colors(styles: list[ChainStyle]) -> list[ChainStyle]:
    """Fill in missing colors from DEFAULT_COLORS."""
    color_iter = iter(DEFAULT_COLORS)
    resolved   = []
    for s in styles:
        c = s.color if s.color is not None else next(color_iter)
        resolved.append(ChainStyle(
            label      = s.label,
            filled     = s.filled,
            fill_alpha = s.fill_alpha,
            line_alpha = s.line_alpha,
            linestyle  = s.linestyle,
            linewidth  = s.linewidth,
            color      = c,
        ))
    return resolved


def _print_constraints(
    samples_list: list[MCSamples],
    styles:       list[ChainStyle],
    param_names:  list[str],
) -> None:
    for style, sample in zip(styles, samples_list):
        print(f"\n── {style.label} ──")
        existing = sample.getParamNames().list()
        for name in param_names:
            if name not in existing:
                continue
            pretty = sample.getParamNames().parWithName(name).label
            latex  = sample.getInlineLatex(name, limit=1)
            display(Math(rf"{pretty}: {latex}"))


# ══════════════════════════════════════════════════════════════════════════════
# Main public function
# ══════════════════════════════════════════════════════════════════════════════

def plot_triangle(
    paths:   list[str],
    params,
    styles:  list[ChainStyle],
    name:    str        = "triangle_plot",
    config:  PlotConfig = None,
    save: bool = True,
    multiple_fiducials: bool = False,
) -> list[MCSamples]:
    """
    Generate a paper-ready GetDist triangle plot.

    Parameters
    ----------
    paths   : One chain file path per chain, matched to `styles`.
    params  : List of (param_name, fiducial_value) tuples.
              If multiple_fiducials=True, this should instead be a list
              of parameter lists, one per chain.
    styles  : One ChainStyle per chain.
    name    : Output filename stem.
    config  : PlotConfig instance.
    save    : Save figure.
    multiple_fiducials : If True, use one fiducial dictionary per chain.

    Returns
    -------
    samples_list : Loaded GetDist MCSamples objects.
    """
    if config is None:
        config = PlotConfig()

    if len(paths) != len(styles):
        raise ValueError(
            f"len(paths)={len(paths)} must equal len(styles)={len(styles)}"
        )

    if multiple_fiducials:
        if len(params) != len(paths):
            raise ValueError(
                f"Expected {len(paths)} parameter lists, got {len(params)}."
            )

        param_names = [p[0] for p in params[0]]
        fiducial_values = [{p[0]: p[1] for p in par} for par in params]
    else:
        param_names = [p[0] for p in params]
        fiducial_values = {p[0]: p[1] for p in params}

    styles = _resolve_colors(styles)

    # ── Load chains ───────────────────────────────────────────────────────────
    samples_list = []
    for path, style in zip(paths, styles):
        s = _load_cosmosis_chain(path, param_names, config.burn_fraction)
        samples_list.append(s)
        print(f"  {style.label}: {len(s.samples):,} samples  ({path})")

    # ── Build GetDist arguments per chain ─────────────────────────────────────
    colors = [s.color for s in styles]
    filled = [s.filled for s in styles]
    alphas = [
        s.fill_alpha if s.filled else s.line_alpha
        for s in styles
    ]
    contour_ls = [s.linestyle for s in styles]
    contour_lws = [s.linewidth for s in styles]
    labels = [s.label for s in styles]

    # ── Render ────────────────────────────────────────────────────────────────
    with mpl.rc_context(_PAPER_RC):
        g = plots.get_subplot_plotter(width_inch=config.figsize)
        g.settings.axes_fontsize = 13
        g.settings.lab_fontsize = 18
        g.settings.legend_fontsize = 15
        g.settings.figure_legend_loc = "upper right"
        g.settings.axis_tick_x_rotation = 45
        g.settings.num_plot_contours = 2

        g.triangle_plot(
            samples_list,
            legend_labels=labels,
            filled=filled,
            colors=colors,
            contour_colors=colors,
            alphas=alphas,
            contour_ls=contour_ls,
            contour_lws=contour_lws,
            fine_bins=1,
            markers=None if multiple_fiducials else fiducial_values,
            param_limits=config.param_limits,
        )

        if multiple_fiducials:
            for marker, style in zip(fiducial_values, styles):
                g.add_param_markers(
                    marker,
                    color=style.color,
                    ls=":",
                    lw=0.8,
                )

        out_path = Path(config.output_dir) / f"{name}.pdf"
        if save:
            plt.savefig(out_path, bbox_inches="tight", dpi=config.dpi)
            print(f"\nSaved to {out_path}")

        plt.show()

    # ── Constraints ───────────────────────────────────────────────────────────
    if config.print_constraints:
        _print_constraints(samples_list, styles, param_names)

    return samples_list

# ══════════════════════════════════════════════════════════════════════════════
# FITS export
# ══════════════════════════════════════════════════════════════════════════════

def chains_to_fits(
    params:        list[tuple[str, float]],
    labels:        list[str],
    output_dir:    str                  = ".",
    paths:         list[str]            = None,
    samples:       list[MCSamples]      = None,
    burn_fraction: float                = 0.3,
    overwrite:     bool                 = True,
) -> list[Path]:
    """
    Save one FITS file per chain.  Each file contains:
      - One column per parameter (float64).
      - Fiducial values stored as header keywords  FID_<PARAM>.
      - Chain label stored as header keyword       LABEL.
      - Source path stored as header keyword       SRC_PATH  (if loaded from file).

    Exactly one of `paths` or `samples` must be provided.

    Parameters
    ----------
    params       : List of (param_name, fiducial_value) — same as plot_triangle.
    labels       : One short label per chain, used to name the output files
                   (e.g. ["full", "lensing", "counts"] → full.fits, …).
    output_dir   : Directory where the FITS files are written.
    paths        : CosmoSIS chain file paths.  Mutually exclusive with `samples`.
    samples      : Already-loaded MCSamples list.  Mutually exclusive with `paths`.
    burn_fraction: Burn-in to remove when loading from `paths` (ignored for `samples`).
    overwrite    : Overwrite existing files silently.

    Returns
    -------
    out_paths : List of Path objects for the written files.
    """
    if (paths is None) == (samples is None):
        raise ValueError("Provide exactly one of `paths` or `samples`, not both.")
    if paths is not None and len(paths) != len(labels):
        raise ValueError(f"len(paths)={len(paths)} must equal len(labels)={len(labels)}")
    if samples is not None and len(samples) != len(labels):
        raise ValueError(f"len(samples)={len(samples)} must equal len(labels)={len(labels)}")

    param_names     = [p[0] for p in params]
    fiducial_values = {p[0]: p[1] for p in params}
    out_dir         = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load from disk if needed
    if paths is not None:
        chain_list = [
            (_load_cosmosis_chain(p, param_names, burn_fraction), p)
            for p in paths
        ]
    else:
        chain_list = [(s, None) for s in samples]

    # Sanitize a param name to a valid FITS column name
    def _fits_col(name: str) -> str:
        return re.sub(r"[^0-9a-zA-Z_]", "_", name).strip("_") or "param"

    out_paths = []
    for (sample, src_path), label in zip(chain_list, labels):
        array = sample.samples  # shape (N, n_params) — already burn-trimmed

        # ── Build astropy Table ───────────────────────────────────────────────
        table = Table()
        for i, name in enumerate(param_names):
            table[_fits_col(name)] = array[:, i].astype(np.float64)

        # ── Header metadata ───────────────────────────────────────────────────
        table.meta["LABEL"]    = label
        table.meta["N_SAMP"]   = len(array)
        table.meta["N_PAR"]    = len(param_names)
        if src_path is not None:
            table.meta["SRC_PATH"] = str(src_path)

        # Fiducial values — FITS keywords are ≤ 8 chars, so we abbreviate
        for name, fid in fiducial_values.items():
            key = f"FID_{_fits_col(name)}"[:8]
            table.meta[key] = fid

        # Raw parameter names stored as comment cards (preserves LaTeX strings)
        for col_name, raw_name in zip([_fits_col(n) for n in param_names], param_names):
            table.meta[f"RAW_{col_name}"[:8]] = raw_name

        # ── Write ─────────────────────────────────────────────────────────────
        out_path = out_dir / f"{label}.fits"
        table.write(str(out_path), overwrite=overwrite)
        print(f"  Saved {len(array):,} samples × {len(param_names)} params → {out_path}")
        out_paths.append(out_path)

    return out_paths


def fits_to_samples(
    fits_paths:  list[str | Path],
    burn_fraction: float = 0.0,
) -> list[MCSamples]:
    """
    Reload chains previously saved by chains_to_fits() back into MCSamples.

    Parameters
    ----------
    fits_paths    : Paths to the FITS files written by chains_to_fits().
    burn_fraction : Additional burn-in to remove on reload (default 0 — already trimmed).

    Returns
    -------
    samples_list : Ready-to-use MCSamples objects.
    """
    samples_list = []
    for fp in fits_paths:
        table  = Table.read(str(fp))
        array  = np.column_stack([table[col].data for col in table.colnames])
        names  = [table.meta.get(f"RAW_{col}"[:8], col) for col in table.colnames]
        label  = table.meta.get("LABEL", Path(fp).stem)

        sample = MCSamples(
            samples  = array,
            names    = names,
            labels   = names,
            settings = GETDIST_SETTINGS,
        )
        if burn_fraction > 0:
            sample.removeBurn(burn_fraction)

        print(f"  Loaded {len(sample.samples):,} samples ({label})  ← {fp}")
        samples_list.append(sample)

    return samples_list


def plot_box_triangle(
    paths_lower:  list[str],
    styles_lower: list[ChainStyle],
    paths_upper:  list[str],
    styles_upper: list[ChainStyle],
    params:       list[tuple[str, float]],
    name:         str = "box_triangle",
    config:       PlotConfig = None,
    save:         bool = True,
) -> tuple[list[MCSamples], list[MCSamples]]:
    if config is None:
        config = PlotConfig()

    param_names     = [p[0] for p in params]
    fiducial_values = {p[0]: p[1] for p in params}

    styles_lower = _resolve_colors(styles_lower)
    styles_upper = _resolve_colors(styles_upper)

    samples_lower = [
        _load_cosmosis_chain(p, param_names, config.burn_fraction)
        for p in paths_lower
    ]
    samples_upper = [
        _load_cosmosis_chain(p, param_names, config.burn_fraction)
        for p in paths_upper
    ]

    with mpl.rc_context(_PAPER_RC):
        g = plots.get_subplot_plotter(width_inch=config.figsize)
        g.settings.axes_fontsize        = 13
        g.settings.lab_fontsize         = 18
        g.settings.legend_fontsize      = 15
        g.settings.axis_tick_x_rotation = 45
        g.settings.num_plot_contours    = 2
        g.settings.figure_legend_loc = "lower left"
        g.triangle_plot(
            samples_lower,
            params=param_names,
            filled=[s.filled for s in styles_lower],
            contour_colors=[s.color for s in styles_lower],
            contour_ls=[s.linestyle for s in styles_lower],
            contour_lws=[s.linewidth for s in styles_lower],
            legend_labels=[s.label for s in styles_lower],
            markers=fiducial_values,
            upper_roots=samples_upper,
            upper_kwargs=dict(
                filled=[s.filled for s in styles_upper],
                contour_colors=[s.color for s in styles_upper],
                contour_ls=[s.linestyle for s in styles_upper],
                contour_lws=[s.linewidth for s in styles_upper],
                show_1d=True,
            ),
            upper_label_right=True,
        )

        # Second legend for the upper-triangle chains (triangle_plot only
        # auto-legends the main `roots`, not upper_roots)
        upper_handles = [
            plt.Line2D([0], [0], color=s.color, lw=s.linewidth, ls=s.linestyle)
            for s in styles_upper
        ]
        g.fig.legend(
            upper_handles, [s.label for s in styles_upper],
            loc="upper center",
            bbox_to_anchor=(0.5, 1.02),   # centered, just above the grid
            ncol=len(styles_upper),       # one row
            frameon=False,
            fontsize=13,
        )
        
        # make room so the top legend doesn't overlap the top row of panels
        g.fig.subplots_adjust(top=0.93)

        out_path = Path(config.output_dir) / f"{name}.pdf"
        Path(config.output_dir).mkdir(parents=True, exist_ok=True)
        if save:
            plt.savefig(out_path, bbox_inches="tight", dpi=config.dpi)
            print(f"\nSaved to {out_path}")
        plt.show()

    if config.print_constraints:
        print("\n=== Lower triangle ===")
        _print_constraints(samples_lower, styles_lower, param_names)
        print("\n=== Upper triangle ===")
        _print_constraints(samples_upper, styles_upper, param_names)

    return samples_lower, samples_upper