## Example of running clustra on a single generated data set
## 
library(clustra)

## Set working directory
playdir = "~/clustra_play"
if(!dir.exists(playdir)) dir.create(playdir)

## Set or get clustra parameters
## TODO document full PL.  cores A list specifying multicore parallelism with
## components e_mc (expectation across k), m_mc (maximization across k),
## bam_nthreads (see bam documentation), blas (OpenBlas). Care should be taken
## that cores are not oversubscribed.
PL = clustra_par(playdir)

## Set seed for reproducibility
## TODO Check if reproducible when running parallel
set.seed(PL$gen_par$seed)
a0 = a = deltime()

## Generate a data set
data = gen_traj_data(PL)
write.table(data, file = "clustra_gen_data.csv")
## plot to check result
iplot = sample(unique(data$id), 8)
sampobs = match(data$id, iplot, nomatch = 0) > 0
ggplot2::ggplot(data[sampobs, ], ggplot2::aes(x = time, y = response)) +
  ggplot2::facet_wrap(~ id) + ggplot2::geom_point()

a = a_fit = deltime(a, paste0("\nData (", paste(dim(data), collapse = ","), ") generated"))
head(data)

## cluster the trajectories
cl = clustra(data, 3, PL, verbose = TRUE, plot = TRUE)
