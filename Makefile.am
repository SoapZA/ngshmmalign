AM_CPPFLAGS = -I$(top_srcdir)/src $(BOOST_CPPFLAGS)

EXTRA_DIST = CMakeLists.txt

dist_doc_DATA = README.md

bin_PROGRAMS = ngshmmalign

ngshmmalign_SOURCES = \
	src/aligner.hpp \
	src/aligner_impl.hpp \
	src/debug.hpp \
	src/dna_array.hpp \
	src/fasta.hpp \
	src/fastq.hpp \
	src/hmmalign.hpp \
	src/hmmalign_impl.hpp \
	src/index_impl.hpp \
	src/reference.hpp \
	src/reference_impl.hpp \
	src/utility_functions.hpp \
	src/main.cpp

ngshmmalign_LDADD = \
	$(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_IOSTREAMS_LIB)

# CHECKS
TESTS = $(check_PROGRAMS)
check_PROGRAMS = \
	testsuite/dna_array_compile_test \
	testsuite/type_caster_compile_test \
	testsuite/reference_compile_test \
	testsuite/hmmalign_forward_compile_test

testsuite_dna_array_compile_test_SOURCES = testsuite/dna_array_compile_test.cpp

testsuite_type_caster_compile_test_SOURCES = testsuite/type_caster_compile_test.cpp

testsuite_reference_compile_test_SOURCES = testsuite/reference_compile_test.cpp
testsuite_reference_compile_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB)

testsuite_hmmalign_forward_compile_test_SOURCES = testsuite/hmmalign_forward_compile_test.cpp
testsuite_hmmalign_forward_compile_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB)
