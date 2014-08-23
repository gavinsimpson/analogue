# Get the version info for later
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

all: check clean

#docs:
#	R -q -e 'library("roxygen2"); roxygenise(".")'

build: #docs
	cd ..;\
	R CMD build analogue

check: build
	cd ..;\
	R CMD check analogue_$(PKGVERS).tar.gz

check-cran: build
	cd ..;\
	R CMD check --as-cran analogue_$(PKGVERS).tar.gz

install: build
	cd ..;\
	R CMD INSTALL analogue_$(PKGVERS).tar.gz

move: check
	cp ../analogue.Rcheck/analogue-Ex.Rout ./tests/Examples/analogue-Ex.Rout.save

clean:
	cd ..;\
	rm -r analogue.Rcheck/
