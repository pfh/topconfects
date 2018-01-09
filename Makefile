
quick : document
	R CMD build --no-build-vignettes .
	R CMD check --no-build-vignettes topconfects_*.tar.gz
	rm topconfects_*.tar.gz

check :
	Rscript -e "devtools::check()"

document :
	Rscript -e "devtools::document()"

site : document
	echo "pkgdown::build_site()" |R --vanilla

publish : 
	scp -r docs/* logarithmic.net:www/topconfects/

