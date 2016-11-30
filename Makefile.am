AM_CPPFLAGS = -I$(top_srcdir)/src $(BOOST_CPPFLAGS)

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
	dna_array_compile_test \
	type_caster_compile_test \
	reference_compile_test \
	hmmalign_forward_compile_test

dna_array_compile_test_SOURCES = testsuite/dna_array_compile_test.cpp
type_caster_compile_test_SOURCES = testsuite/type_caster_compile_test.cpp

reference_compile_test_SOURCES = testsuite/reference_compile_test.cpp
reference_compile_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB)

hmmalign_forward_compile_test_SOURCES = testsuite/hmmalign_forward_compile_test.cpp
hmmalign_forward_compile_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB)
