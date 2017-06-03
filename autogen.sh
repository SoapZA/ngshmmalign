#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	rm -rf .deps/
	rm -rf .Tpo
	rm -rf *.dSYM
	rm -rf ngshmmalign
	rm -rf ngshmmalign-*.tar.bz2

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
	rm -rf *_compile_test
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
