dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[ngshmmalign],
	[0.1.1],
	[david.seifert@bsse.ethz.ch],
	[ngshmmalign],
	[https://github.com/cbg-ethz/ngshmmalign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
AC_MSG_CHECKING([whether configure should try to set CXXFLAGS/CPPFLAGS/LDFLAGS])
AS_IF([test "x${CXXFLAGS+set}" = "xset" || test "x${CPPFLAGS+set}" = "xset" || test "x${LDFLAGS+set}" = "xset"], [
	enable_flags_setting=no
], [
	enable_flags_setting=yes
])
AC_MSG_RESULT([${enable_flags_setting}])

AX_CHECK_ENABLE_DEBUG
AS_IF([test "x${enable_debug}" = "xno"], [
	cxxflags_test="-ffast-math -pedantic -O3"
	cppflags_test=""

	AC_CANONICAL_HOST
	AS_CASE([${host_os}],
		[darwin*],
			[ldflags_test="-Wl,-dead_strip_dylibs"],
		[linux*],
			[ldflags_test="-Wl,-O1 -Wl,--as-needed"]
	)
])



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AX_COMPILER_VENDOR
AX_COMPILER_VERSION
AC_PROG_SED
AC_PROG_CXX



dnl ================
dnl Check for OpenMP
dnl ================
AC_OPENMP
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF([test "x${enable_flags_setting}" = "xyes" && test "x${enable_debug}" = "xno"], [
	AX_APPEND_COMPILE_FLAGS([${cxxflags_test} -Wall -Wno-redeclared-class-member -Wno-unused-local-typedef], [CXXFLAGS])
	CXXFLAGS=$( echo ${CXXFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

	AX_APPEND_COMPILE_FLAGS([${cppflags_test}], [CPPFLAGS])
	CPPFLAGS=$( echo ${CPPFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

	AX_APPEND_LINK_FLAGS([${ldflags_test}], [LDFLAGS])
	LDFLAGS=$( echo ${LDFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )
])



dnl ===============
dnl Check for C++11
dnl ===============
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])



dnl ===============
dnl Check for Boost
dnl ===============
AX_BOOST_BASE([1.59], [], [AC_MSG_ERROR([could not find a suitable boost installation.])])
AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_FILESYSTEM
AX_BOOST_IOSTREAMS



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.14 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	${PACKAGE_NAME} ${VERSION}
	Using ${ax_cv_cxx_compiler_vendor}, ${ax_cv_cxx_compiler_version}

	CXX:        ${CXX}
	CXXFLAGS:   ${CXXFLAGS}
	CPPFLAGS:   ${CPPFLAGS}
	LDFLAGS:    ${LDFLAGS}
])
