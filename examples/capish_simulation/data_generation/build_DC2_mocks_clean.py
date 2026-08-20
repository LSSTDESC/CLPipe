import numpy as np
import jaxlib, jax
import jax.numpy as jnp
import equinox, optax
import matplotlib.pyplot as plt
import sys
import os
import pyccl as ccl
import classy
import PySSC
import pyccl
from astropy.table import Table
from pathlib import Path
# Add project root to path
sys.path.insert(0, "/sps/lsst/users/ebarroso/capish/")
from modules.simulation import UniverseSimulator
import configparser
import clmm
from clmm.dataops import compute_tangential_and_cross_components, make_radial_profile, make_bins
from clmm.galaxycluster import GalaxyCluster
import clmm.utils as u
from clmm import Cosmology
from clmm.support import mock_data as mock

INI_PATH = "capish_config_DC2like_clean.ini"
OUTPUT_DIR = "../../../data/capish_simulation/data_generation/mocks_seeds/"

code, seed1, seed2 = sys.argv[0], int(sys.argv[1]), int(sys.argv[2])

os.makedirs(OUTPUT_DIR, exist_ok=True)


def ccl_cosmo(config_new):
    Omega_m = float(config_new['parameters']['Omega_m'])
    Omega_b = float(config_new['parameters']['Omega_b'])
    sigma8 = float(config_new['parameters']['sigma8'])
    h = float(config_new['parameters']['h'])
    ns = float(config_new['parameters']['ns'])
    w0 = float(config_new['parameters']['w0'])
    wa = float(config_new['parameters']['wa'])
    cosmo_ccl_object = ccl.Cosmology(
        Omega_c=Omega_m - Omega_b, Omega_b=Omega_b, h=h, sigma8=sigma8, n_s=ns,
        w0=w0, wa=wa, transfer_function='boltzmann_class', matter_power_spectrum='linear',
    )
    cosmo_clmm_object = Cosmology(H0=100 * h, Omega_dm0=Omega_m - Omega_b, Omega_b0=Omega_b, Omega_k0=0.0)
    cosmo_clmm_object.be_cosmo = cosmo_ccl_object
    return cosmo_ccl_object, cosmo_clmm_object


default_config_capish = configparser.ConfigParser()
default_config_capish.read(INI_PATH)

Omega_m = float(default_config_capish['parameters']['Omega_m'])
sigma8 = float(default_config_capish['parameters']['sigma8'])

shape_noise = None
ngal_density = float(default_config_capish['cluster_catalogue']['ngal_arcmin2'])
concentration_fixed_value = float(default_config_capish['cluster_catalogue']['concentration_fixed_value'])

r_min = float(default_config_capish['cluster_catalogue']['DeltaSigma_Rmin'])
r_max = float(default_config_capish['cluster_catalogue']['DeltaSigma_Rmax'])
n_radius_bins = int(default_config_capish['cluster_catalogue']['n_radius_bins'])
radius_bin_method = default_config_capish['cluster_catalogue']['radius_bin_method']

richness_edges = [float(x) for x in default_config_capish['summary_statistics']['richness_edges'].split(',')]
redshift_edges = [float(x) for x in default_config_capish['summary_statistics']['redshift_edges'].split(',')]
richness_min, richness_max = min(richness_edges), max(richness_edges)
z_min, z_max = min(redshift_edges), max(redshift_edges)

field_size_mpc = 7.0

cosmo_ccl_object, cosmo_clmm_object = ccl_cosmo(default_config_capish)


def generate_profile(log10m, z, c, cosmo, new_bins):
    mdelta = 10 ** log10m
    noisy_data_z = mock.generate_galaxy_catalog(
        mdelta, z, c,
        cosmo, "chang13", zsrc_min=z + 0.2,
        massdef="critical", shapenoise=shape_noise, photoz_sigma_unscaled=None,
        ngal_density=ngal_density, cluster_ra=0, cluster_dec=0,
        field_size=field_size_mpc,
    )
    cl = GalaxyCluster('id', 0, 0, z, noisy_data_z)
    cl.compute_tangential_and_cross_components(add=True, cosmo=cosmo, is_deltasigma=True)
    r_proj_i = u.convert_units(
        cl.galcat['theta'], "radians", "Mpc", redshift=z, cosmo=cosmo
    )
    kappa_i = clmm.compute_convergence(
        np.asarray(r_proj_i), mdelta, c, z, np.asarray(cl.galcat['z']), cosmo,
        delta_mdef=200, halo_profile_model="nfw", massdef="critical",
    )
    cl.galcat['et'] = cl.galcat['et'] * (1.0 - kappa_i)
    cl.galcat['ex'] = cl.galcat['ex'] * (1.0 - kappa_i)

    cl.compute_galaxy_weights(
        use_pdz=False, use_shape_noise=False, use_shape_error=False,
        is_deltasigma=True, cosmo=cosmo,
    )
    new_profiles = cl.make_radial_profile(
        "Mpc", bins=new_bins, cosmo=cosmo, use_weights=True, weights_in="w_ls"
    )
    return new_profiles['radius'], new_profiles['gt'], new_profiles['W_l']


sim = UniverseSimulator(
    default_config_path=None, default_config=default_config_capish,
    variable_params_names=['Omega_m', 'sigma8'],
)

np.random.seed(seed1)
log10m_halo, z_true, richness, log10mWL, z_obs = sim.run_simulation_halo_and_cluster_catalogue(
    [Omega_m, sigma8]
)

mask_cluster_catalog = (
    (z_obs >= z_min) * (z_obs <= z_max)
    * (richness >= richness_min) * (richness <= richness_max)
)

new_bins = make_bins(r_min, r_max, nbins=n_radius_bins, method=radius_bin_method)

radius_list, DS_list, Wl_list = [], [], []
indexes = np.arange(len(log10m_halo))[mask_cluster_catalog]

for index in indexes:
    np.random.seed(seed2 + index)
    r, DS, Wl = generate_profile(
        log10m_halo[index], z_true[index], concentration_fixed_value, cosmo_clmm_object, new_bins
    )
    radius_list.append(r)
    DS_list.append(DS)
    Wl_list.append(Wl)

Summary_table = Table()
Summary_table['binned_radius_Mpc'] = np.array(radius_list)
Summary_table['binned_profiles'] = np.array(DS_list)
Summary_table['binned_weights'] = np.array(Wl_list)
Summary_table['cluster_redshift'] = z_obs[indexes]
Summary_table['cluster_richness'] = richness[indexes]

out_path = os.path.join(OUTPUT_DIR, f"mock_DC2like_data_shape_seed{seed1}_{seed2}.fits")
Summary_table.write(out_path, overwrite=True)
print(f"wrote {len(indexes)} clusters to {out_path}")
