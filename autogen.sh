#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	rm -rf *.dSYM
	rm -rf ngshmmalign
	rm -rf ngshmmalign-0.1.tar.bz2

	# Autotools
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf compile
	rm -rf config.*
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf stamp-h1
	rm -rf test-driver
	
	# source
	rm -rf src/*.o
	rm -rf src/.deps
	rm -rf src/.dirstamp
	
	# testsuite
	rm -rf *.log
	rm -rf *.trs
	rm -rf dna_array_compile_test
	rm -rf hmmalign_forward_compile_test
	rm -rf type_caster_compile_test
	rm -rf parameter_pack_compile_test
	rm -rf testsuite/*.o
	rm -rf testsuite/.deps
	rm -rf testsuite/.dirstamp
	
	# OS X cruft
	find . -name '.DS_Store' -type f -delete
}

if [[ "$1" == "--clean" ]]
then
	clean_files
	exit
fi

echo "Bootstrapping Autotools"
autoreconf -vif

if [[ "$1" == "--test" ]]
then
	echo "${DISTCHECK_CONFIGURE_FLAGS}"
	./configure ${DISTCHECK_CONFIGURE_FLAGS}
	DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}" make distcheck
	clean_files
fi
