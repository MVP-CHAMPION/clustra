# A package for clustering longitudinal trajectories. 

Clusters medical trajectories (time series) aligned by an intervention at time 0. Clustering proceeds by an EM algorithm that iterates switching between fitting a bspline to combined responses within each cluster (M-step) and reassigning cluster membership based on nearest fitted bspline (E-step). Initial cluster assignments are random. The fitting is done with the *mgcv* package function *bam*, which scales well to very large data sets.

See the vignette for detailed use examples.