import os
import pyccl as ccl
import sacc
from firecrown.likelihood.gaussian import ConstGaussian
from firecrown.likelihood.binned_cluster_number_counts import BinnedClusterNumberCounts
from firecrown.likelihood.likelihood import Likelihood, NamedParameters
from firecrown.modeling_tools import ModelingTools
from firecrown.models.cluster.abundance import ClusterAbundance
from firecrown.models.cluster.properties import ClusterProperty
from firecrown.models.cluster.recipes.murata_binned_spec_z import MurataBinnedSpecZRecipe

def get_cluster_abundance() -> ClusterAbundance:
    '''Creates and returns a ClusterAbundance object.''' 
    hmf = ccl.halos.MassFuncDespali16(mass_def='200c')  # Using despali16 from the config
    min_mass, max_mass = 13.0, 16.0
    min_z, max_z = 0.2, 0.8
    cluster_abundance = ClusterAbundance((min_mass, max_mass), (min_z, max_z), hmf)

    return cluster_abundance

def build_likelihood(build_parameters: NamedParameters) -> tuple[Likelihood, ModelingTools]:
    '''Builds the likelihood for Firecrown.''' 
    # Pull params for the likelihood from build_parameters
    average_on = ClusterProperty.NONE
    if build_parameters.get_bool('use_cluster_counts', True):
        average_on |= ClusterProperty.COUNTS
    if build_parameters.get_bool('use_mean_log_mass', True):
        average_on |= ClusterProperty.MASS
    if build_parameters.get_bool('use_mean_deltasigma', True):
        average_on |= ClusterProperty.DELTASIGMA

    recipe_counts = MurataBinnedSpecZRecipe()
    recipe_counts.mass_distribution.pivot_mass = 32.79792227308491
    recipe_counts.mass_distribution.pivot_redshift = 0.15
    recipe_counts.mass_distribution.log1p_pivot_redshift = 0.13976194237515868
    survey_name = 'SDSSCluster_redshift_richness'
    likelihood = ConstGaussian(
        [BinnedClusterNumberCounts(average_on, survey_name, recipe_counts)]
    )

    sacc_path = 'clusters_SDSS_sacc_file_cov.sacc'
    sacc_data = sacc.Sacc.load_fits(sacc_path)
    likelihood.read(sacc_data)

    cluster_abundance = get_cluster_abundance()
    modeling_tools = ModelingTools(cluster_abundance=cluster_abundance)

    return likelihood, modeling_tools
