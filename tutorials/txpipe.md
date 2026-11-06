# Measurements with TXPipe

CLPipe does not measure the cluster observables itself. The number counts and
the stacked $\Delta\Sigma$ profiles are measured by
[TXPipe](https://github.com/LSSTDESC/TXPipe), the DESC measurement pipeline,
which has a set of stages for clusters. A CLPipe run calls them through the
`TXPipe.yml` file of its run directory, in the `txpipe_clp` environment, and
the next steps start from the SACC file they write.

The code and the documentation of these stages are in TXPipe
(https://txpipe.readthedocs.io). This page only says what CLPipe uses and
where to look.

## Version

The measurements of the CLPipe paper were made with the TXPipe tag
[`clpipe-paper-v1`](https://github.com/LSSTDESC/TXPipe/tree/clpipe-paper-v1).
Use it to reproduce them.

## The stages

The cluster stages are in
[`txpipe/extensions/cluster_counts`](https://github.com/LSSTDESC/TXPipe/tree/clpipe-paper-v1/txpipe/extensions/cluster_counts).
They run one after the other, each reading the output of the previous one.

| Stage | Output | What it does |
|---|---|---|
| `CLClusterBinningRedshiftRichness` | `cluster_catalog_tomography` | Splits the cluster catalog in bins of redshift and richness |
| `CLClusterShearCatalogs` | `cluster_shear_catalogs` | Finds the background sources around each cluster and computes their weights |
| `CLClusterEnsembleProfiles` | `cluster_profiles` | Stacks the $\Delta\Sigma$ profiles in each bin with CLMM and estimates their covariance |
| `CLClusterSACC` | `cluster_sacc_catalog` | Writes the counts, the profiles and the profile covariance in a SACC file |

Before them, `TXSourceSelectorMetadetect` selects the source galaxies and
calibrates their shear, using the tomographic classifier trained by
`TXSourceTomography`. The photometric redshifts of the sources come from the
RAIL stages `BPZliteInformer` and `BPZliteEstimator`.

## Configuration

Each stage has a section in the configuration file of the run, next to the
sections of the CLPipe stages. In the cosmoDC2 redMaPPer analysis
(`examples/cosmodc2_redmapper/baseline/cosmodc2_redmapper_full_analysis`), for
instance

```yaml
CLClusterBinningRedshiftRichness:
    zedge: [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    richedge: [20.0, 35.0, 70.0, 100.0, 200.0]

CLClusterEnsembleProfiles:
    r_min: 1.1      # Mpc
    r_max: 3.68
    nbins: 7
    cov_type: bootstrap_covariance

CLClusterSACC:
    survey_name: cosmodc2_redmapper
    area: 440.0     # deg2
    cov_component: tan_bs   # bootstrap covariance of the tangential profiles
```

## What CLPipe reads

`cluster_sacc_catalog` is the `clusters_sacc_file` input of `CLPCovariance`.
It holds the counts in each redshift and richness bin, the profiles in each
radial bin, and the covariance of the profiles. `CLPCovariance` adds the
covariance of the counts (see
[covariance_tjpcov.ipynb](covariance_tjpcov.ipynb)).

To try the stages on a small sample, the `examples/cosmoDC2-1deg2` run uses
the 1 deg$^2$ TXPipe example data, and its README explains how to get it.
