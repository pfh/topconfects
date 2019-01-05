
quick : document
	R CMD build --no-build-vignettes .
	R CMD check --no-build-vignettes topconfectswald_*.tar.gz
	rm topconfectswald_*.tar.gz

check :
	Rscript -e "devtools::check()"

rmd_tests :
	Rscript -e "rmarkdown::render('test_rmds/test_limma.Rmd')"
	Rscript -e "rmarkdown::render('test_rmds/test_edgeR.Rmd')"

document :
	Rscript -e "devtools::document()"

install : document
	Rscript -e "devtools::install()"

#site : install
#	echo "pkgdown::build_site()" |R --vanilla

#publish : 
#	scp -r docs/* logarithmic.net:www/topconfects/

