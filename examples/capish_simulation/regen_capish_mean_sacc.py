"""Standalone reproduction of analysis_summary.ipynb cells 2 + 12: builds the
mean data vector and its realization-scatter covariance from the 1000 capish
mock seed realizations, and saves mock_DC2like_cluster_sacc_mean.sacc.

Companion to regen_crow_sacc.py (which saves the crow theory prediction as
its own sacc file). Run this whenever the underlying mock seeds change.
"""
import glob
import numpy as np
from astropy.table import Table
import sacc
from scipy.linalg import block_diag

# ---- load 1000 seed realizations, build mean + covariance ----
seed_files = sorted(glob.glob("../../data/capish_simulation/data_generation/mocks_seeds/mock_DC2like_data_shape_seed*.fits"))
print(f"Found {len(seed_files)} seed realizations")

def binning(corner):
    return [[corner[i], corner[i + 1]] for i in range(len(corner) - 1)]

z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
Z_bin = np.array(binning(z_corner))
rich_corner = np.array([20, 35, 70, 100, 200])
Obs_bin = np.array(binning(rich_corner))

n_seeds = len(seed_files)
n_z = len(Z_bin)
n_rich = len(Obs_bin)

counts_all = np.zeros((n_seeds, n_z, n_rich))
lensing_all = None
radius_ref = None

for s_idx, fpath in enumerate(seed_files):
    lensing_profiles = Table.read(fpath)

    if radius_ref is None:
        n_radius = lensing_profiles["binned_radius_Mpc"].shape[1]
        radius_ref = np.mean(lensing_profiles["binned_radius_Mpc"], axis=0)
        lensing_all = np.zeros((n_seeds, n_z, n_rich, n_radius))

    for i, z_bin in enumerate(Z_bin):
        for j, rich_bin in enumerate(Obs_bin):
            mask = (lensing_profiles["cluster_redshift"] > z_bin[0]) * (lensing_profiles["cluster_redshift"] <= z_bin[1])
            mask *= (lensing_profiles["cluster_richness"] > rich_bin[0]) * (lensing_profiles["cluster_richness"] <= rich_bin[1])

            counts_all[s_idx, i, j] = np.sum(mask)

            if np.sum(mask) > 0:
                lensing_profiles_bin = lensing_profiles[mask]
                lensing_all[s_idx, i, j] = np.average(
                    lensing_profiles_bin["binned_profiles"],
                    weights=lensing_profiles_bin["binned_weights"],
                    axis=0,
                )

counts_mean = np.mean(counts_all, axis=0)
lensing_mean = np.mean(lensing_all, axis=0)

counts_flat = counts_all.reshape(n_seeds, n_z * n_rich)
counts_cov = np.cov(counts_flat, rowvar=False)

lensing_cov = np.zeros((n_z, n_rich, n_radius, n_radius))
for i in range(n_z):
    for j in range(n_rich):
        vecs = lensing_all[:, i, j, :]
        lensing_cov[i, j] = np.cov(vecs, rowvar=False)

print("counts_mean shape:", counts_mean.shape)
print("lensing_mean shape:", lensing_mean.shape)

# ---- save mean data vector + realization-scatter covariance to sacc ----
survey_name = "capish_dc2"
s = sacc.Sacc()
area = 436.6

for i, zbin in enumerate(Z_bin):
    s.add_tracer("bin_z", f"bin_z_{i}", zbin[0], zbin[1])

for j, rbin in enumerate(Obs_bin):
    s.add_tracer("bin_richness", f"bin_rich_{j}", np.log10(rbin[0]), np.log10(rbin[1]))

for k, radius in enumerate(radius_ref):
    s.add_tracer(
        "bin_radius", f"radius_{k}", radius, radius, radius,
        metadata={"bin_units": "Mpc", "center": radius},
    )

s.add_tracer("survey", survey_name, area)

for i in range(len(Z_bin)):
    for j in range(len(Obs_bin)):
        s.add_data_point(
            sacc.standard_types.cluster_counts,
            (survey_name, f"bin_rich_{j}", f"bin_z_{i}"),
            counts_mean[i, j],
        )

for i in range(len(Z_bin)):
    for j in range(len(Obs_bin)):
        for k in range(n_radius):
            s.add_data_point(
                sacc.standard_types.cluster_delta_sigma,
                (survey_name, f"bin_rich_{j}", f"bin_z_{i}", f"radius_{k}"),
                lensing_mean[i, j, k],
            )

counts_cov_matrix = counts_cov.reshape(len(Z_bin) * len(Obs_bin), len(Z_bin) * len(Obs_bin))
ds_blocks = []
for i in range(len(Z_bin)):
    for j in range(len(Obs_bin)):
        ds_blocks.append(lensing_cov[i, j])

full_cov = block_diag(counts_cov_matrix, *ds_blocks)
s.add_covariance(full_cov)
s.to_canonical_order()

s.save_fits("../../data/capish_simulation/mock_DC2like_cluster_sacc_mean.sacc", overwrite=True)

print("Saved mock_DC2like_cluster_sacc_mean.sacc with area =", area)
print("Number of data points:", len(s.data))
print("Covariance shape:", full_cov.shape)
