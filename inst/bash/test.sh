R --no-save --quiet -e 'devtools::document()'; R CMD build .; R CMD check ShadowCAT_0.1.tar.gz --no-manual --no-build-vignettes
