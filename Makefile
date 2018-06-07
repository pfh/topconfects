
quick : document
	R CMD build --no-build-vignettes .
	R CMD check --no-build-vignettes topconfects_*.tar.gz
	rm topconfectsql_*.tar.gz

check :
	Rscript -e "devtools::check()"

rmd_tests :
	Rscript -e "rmarkdown::render('test_rmds/test_limma.Rmd')"
	Rscript -e "rmarkdown::render('test_rmds/test_edgeR.Rmd')"
	Rscript -e "rmarkdown::render('test_rmds/test_shift_effect.Rmd')"

document :
	Rscript -e "devtools::document()"

#site : document
#	echo "pkgdown::build_site()" |R --vanilla

#publish : 
#	scp -r docs/* logarithmic.net:www/topconfects/

# Documentation file building:
# devtools:::build_manual(path=".")
# devtools:::build_vignettes()
