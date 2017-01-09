
quick :
	Rscript -e "devtools::check(build_args=c('--no-build-vignettes'))"

check :
	Rscript -e "devtools::check()"

document :
	Rscript -e "devtools::document()"

site : document
	Rscript -e "pkgdown::build_site()"

publish : 
	scp -r docs/* logarithmic.net:www/topconfects/

