******************
Rapid overview
******************

clpipe is an automated framework developed within the Legacy
Survey of Space and Time Dark Energy Science Collaboration (LSST DESC) for galaxy cluster cosmology,
in preparation for the first Rubin LSST Data releases. This cluster cosmology pipeline
integrates DESC-developed tools spanning galaxy/cluster catalog ingestion, cluster number count and
stacked weak lensing measurements, modeling, and parameter inference through posterior sampling,
enabling reproducible end-to-end cluster cosmology analyses. We test/validate the pipeline using mock
datasets and the LSST DESC cosmoDC2 simulation, performing a joint analysis of redMaPPer cluster
abundance and stacked weak lensing profiles measured in bins of redshift and richness. We constrain
some cosmological parameters along with six cluster mass–observable relation parameters. Besides
our baseline recovering the fiducial cosmoDC2 values within 2σ, we also showcase the flexibility of the
overall pipeline by assessing the impact of fiducial cosmology, alternative modeling prescriptions, and
cluster finder selection function on the parameter constraints.

.. figure:: diagram.png
   :scale: 100 %
   :alt: CLPipe and its dependencies diagram
   :align: center