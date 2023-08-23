library(clustra)
start_knit = deltime()
library(data.table); library(magrittr);library(ggplot2); library(ggpubr)
library(haven); library(parallel);  library(dplyr)
library(mgcv)
#library(clusteval)
mc = 5
# If running on a unix or a Mac platform, change to number of cores (up to 2x cores on Intel chips with hyperthreading)
# mc = detectCores()
print(paste("Number of cores being used: mc parameter =", mc))
data(bp) # get the package bp data set
data <- bp
head(data)
plot_path <- "./" # output path for plots

set.seed(12345)
plot_sample(data, layout = c(3, 3), group = "true_group")

plot_smooths(data, fits = NULL, max.data = 20000, group = "true_group")

set.seed(12345)
cl = clustra(data, k = 5, maxdf = 30, conv = c(20, 0.5), mccores = mc, verbose = TRUE)

# What happens if we cut the iterations at 5 maximum?
set.seed(12345)
cl.v2 = clustra(data, k = 5, maxdf = 30, conv = c(5, 0.5), mccores = mc, verbose = FALSE)

set.seed(12345)
png(filename = paste0(plot_path, "R_iter%02d_plot.png"))
  cla = clustra(data, k = 5, maxdf = 30, conv = c(10, 0.5), mccores = mc, verbose = 2, ylim = c(110, 170))
dev.off()
knitr::include_graphics(paste0(plot_path, "R_iter01_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter02_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter03_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter04_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter05_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter06_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter07_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter08_plot.png"))
knitr::include_graphics(paste0(plot_path, "R_iter08_plot.png"))

plot_smooths(data, fits = cl$tps, max.data = 30000)
m = matrix(c(0.1, .95, 0.1, .95), nrow = 1)
sc = split.screen(m)
screen(1)
  sc = split.screen(c(1, 2))
  screen(2)
    plot_smooths(data, fits = cl$tps, max.data = 0, ylim = c(110, 180))
  screen(3)
    plot_smooths(data, fits = cl.v2$tps, max.data = 0, ylim = c(110, 180))
close.screen(all.screens = TRUE)

png("R_cl5_plot.png")
    plot_smooths(data, fits = cl$tps, max.data = 0, ylim = c(110, 180))
dev.off()
knitr::include_graphics(paste0(plot_path, "R_cl5_plot.png"))

MixSim::RandIndex(cl$data_group, data[, true_group])$AR # compare the 20-max-iterations to true assignemnt
MixSim::RandIndex(cl.v2$data_group, data[, true_group])$AR # compare the 10-max-iterations to true assignment
MixSim::RandIndex(cl$data_group, cl.v2$data_group)$AR # compare the 10-max-iterations and the 20-max-iterations assignments
#Save the results to the data
data$clus5 <- cl$data_group
data$clus5.v2 <- cl.v2$data_group

set.seed(12345)
t = deltime()
cl10 = clustra(data, k = 10, maxdf = 30, conv = c(50, 0.5), mccores = mc, verbose = TRUE)
deltime(t, paste("run time on", mc, "cores "), units = TRUE, nl = TRUE)
cat(paste("Number of iterations completed:", cl10$iterations))
cat(paste("Number of individuals changing groups in last iteration:", cl10$changes))

png("R_cl10_plot.png")
  plot_smooths(data, cl10$tps, max.data = 0, ylim = c(100, 180))
dev.off()
knitr::include_graphics("R_cl10_plot.png")

MixSim::RandIndex(cl10$data_group, data[, true_group])$AR
data$clus10 <- cl10$data_group # Save the results to the data

set.seed(12345)
t = deltime()
cl2 = clustra(data, k = 2, maxdf = 30, conv = c(20,0.5), mccores=mc, verbose = FALSE)
deltime(t, paste("Seconds run time on", mc, "cores "), units = TRUE, nl = TRUE) # (1.2 minutes on 8 cores)
cat(paste("Number of iterations completed:",cl2$iterations))
cat(paste("Number of individuals changing groups in last iteration:",cl2$changes))
png(paste0(plot_path, "R_cl2_plot.png"))
  plot_smooths(data, cl2$tps, max.data = 0)
dev.off()
knitr::include_graphics("R_cl2_plot.png")
MixSim::RandIndex(cl2$data_group, data[, true_group])$AR
data$clus2 <- cl2$data_group  # Save the results to the data

MixSim::RandIndex(cl$data_group,cl2$data_group)$AR # between 5 and 2 clusters
MixSim::RandIndex(cl$data_group,cl10$data_group)$AR # between 5 and 10 clusters
MixSim::RandIndex(cl10$data_group,cl2$data_group)$AR # between 10 and 2 clusters

t = deltime()
set.seed(12345)
sil = clustra_sil(cl, mccores = mc, conv=c(20,0.5), verbose = FALSE)
lapply(sil, plot_silhouette)
set.seed(12345)
sil2 = clustra_sil(cl2, mccores = mc, conv=c(20,0.5), verbose = FALSE)
lapply(sil2, plot_silhouette)
set.seed(12345)
sil10 = clustra_sil(cl10, mccores = mc, conv=c(20,0.5), verbose = FALSE)
lapply(sil10, plot_silhouette)
t = deltime(t, paste("Silhouettes on", mc, "cores "), units = TRUE, nl = TRUE)

# Example running from the data step
set.seed(12345)
sil = clustra_sil(data, k = c(3,4,20), mccores = mc, conv = c(50,0.5), 
                  verbose = FALSE)
lapply(sil, plot_silhouette)
t = deltime(t, paste("Silhouettes and clustra on", mc, "cores "), units = TRUE, 
            nl = TRUE)

mc = 20
set.seed(12345)
t = deltime()
ran = clustra_rand(data, k = seq(2, 10, 1), starts = "random", mccores = mc, 
                   replicates = 10, conv=c(40,0.5), verbose = TRUE, save = TRUE)
t = deltime(t, paste("Seconds clustra_rand on", mc, "cores "), units = TRUE, 
            nl = TRUE) # 4c-i7 mac 5 cores 3947 seconds
rand_plot(ran, name = "R_randplot.pdf") # save pdf version
rand_plot(ran) # render png version

deltime(start_knit, "clustra vignette run time ", units = TRUE, nl = TRUE)

