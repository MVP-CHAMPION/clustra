## Test environments
* local x86_64-apple-darwin17.0 (64-bit), OS X 10.15.7, R 4.1.0
* win-builder (release and devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'George Ostrouchov <ostrouchovg@ornl.gov>'

New submission

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission
This is a resubmission. In this version, the following are my changes and explanations:

* References: At this time, we don't have a reference that describes our methods, although we are preparing a manuscript that will be available for future releases.  
* All .Rd files now have a `\value` field.  
* `on.exit()` now added to prevent user `par()` changes in `R/evaluate.R`.
* Oak Ridge National Laboratory added as `"cph"` in `Author@R`

Many thanks for your help that now results in a far better documented package!
