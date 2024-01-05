# clustra v0.2.1
* Fixed `difftime` bug for old release R
* Updated maintainer email

# clustra v0.2.0
* Added a 10,000 id data set `bp10k`
* Added a second vignette to reproduce associated paper graphics and more
* Reduced vignette dependence to only base graphics
* New function `plot_sample` built on base graphics
* New function `plot_smooths` built on base graphics
* New function `plot_silhouette` built on base graphics
* Expanded capability of `gen_traj_data` beyond 3 clusters with `type` and `intercept` parameters
* Modified parameters in `clustra`
* Added `starts = "distant"` option in function `start_groups`
* Added parameter `starts` in `clustra_sil` and `clustra_rand`
* New internal functions `kchoose` and `ic_fun` (not exported)

# clustra v0.1.6
* Original ids are now included in `clustra` output
* Added the ability to reuse previous clustra runs by `clustra_sil` for silhouette plots
* `iter` parameter changed to `conv` = c(iter, minchange)

# clustra v0.1.5
Initial release on CRAN