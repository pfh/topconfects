
# Obtain check.Renviron with:
# wget https://raw.githubusercontent.com/Bioconductor/packagebuilder/master/check.Renviron

quick : document
	rm -f topconfects_*.tar.gz
	R CMD build --no-build-vignettes .
	R_CHECK_ENVIRON=check.Renviron R CMD check --no-build-vignettes topconfects_*.tar.gz
	rm topconfects_*.tar.gz

check : document
	rm topconfects_*.tar.gz
	R CMD build .
	R_CHECK_ENVIRON=check.Renviron R CMD check topconfects_*.tar.gz
	rm topconfects_*.tar.gz

bioccheck : document
	#R CMD BiocCheck .
	Rscript -e "BiocCheck::BiocCheck()"

document :
	Rscript -e "devtools::document()"

test :
	Rscript -e "devtools::test()"

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes()"

install : document
	Rscript -e "devtools::install()"

