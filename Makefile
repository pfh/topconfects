
# Obtain check.Renviron with:
# wget https://raw.githubusercontent.com/Bioconductor/packagebuilder/master/check.Renviron

quick : document
	R CMD build --no-build-vignettes .
	R_CHECK_ENVIRON=check.Renviron R CMD check --no-build-vignettes topconfects_*.tar.gz
	rm topconfects_*.tar.gz

check : document
	R CMD build .
	R_CHECK_ENVIRON=check.Renviron R CMD check topconfects_*.tar.gz
	rm topconfects_*.tar.gz

#check :
#	Rscript -e "devtools::check()"

bioccheck : document
	R CMD BiocCheck .

document :
	Rscript -e "devtools::document()"

test :
	Rscript -e "devtools::test()"

install : document
	Rscript -e "devtools::install()"

site : install
	echo "pkgdown::build_site()" |R --vanilla

publish : 
	scp -r docs/* logarithmic.net:www/topconfects/

