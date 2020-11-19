## Example from mgcv::bam to illustrate bam-veclib hang on a mac
## short test
sessionInfo()
set.seed(3)
dat1 <- mgcv::gamSim(1, n=25000, dist="normal", scale=20)
dat2 <- mgcv::gamSim(1, n=25000, dist="normal", scale=20)
bs <- "cr";k <- 12

b1 <- mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat1)
b2 <- mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat2)

dat <- list(dat1, dat2)
b <- lapply(dat, function(dat, bs, k) 
  mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat),
  bs=bs, k=k)
all.equal(b1$coefficients, b[[1]]$coefficients)
all.equal(b2$coefficients, b[[2]]$coefficients)

b <- parallel::mclapply(dat, function(dat, bs, k) 
  mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat),
  bs=bs, k=k, mc.cores = 2)

all.equal(b1$coefficients, b[[1]]$coefficients)
all.equal(b2$coefficients, b[[2]]$coefficients)
## works on Rhea!!
## could it be vecLib??
## try OpenBLAS on Rhea! -> works but no multicore benefit
## try setting single core vecLib on the mac!
