
check :
	Rscript -e "devtools::check(quiet=TRUE)"

document :
	Rscript -e "devtools::document()"

site : document
	Rscript -e "pkgdown::build_site()"

publish : site
	scp -r docs/* logarithmic.net:www/topconfects/

