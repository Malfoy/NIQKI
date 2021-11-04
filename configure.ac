# Must init the autoconf setup
# The first parameter is project name
# second is version number
# third is bug report address
AC_INIT([niqki], [1.0])

# Safety checks in case user overwritten --srcdir
AC_CONFIG_SRCDIR([src/niqki.cpp])

# Store the auxiliary build tools (e.g., install-sh, config.sub, config.guess)
# in this dir (build-aux)
AC_CONFIG_AUX_DIR([config])
dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config/config.h])

# Init automake, and specify this program use relaxed structures.
# i.e. this program doesn't follow the gnu coding standards, and doesn't have
# ChangeLog, COPYING, AUTHORS, INSTALL, README etc. files.
AM_INIT_AUTOMAKE([-Wall -Werror foreign]) 

AC_CHECK_LIB(zlib, zlib)

dnl Save existing environment variables
ORIG_LDADD=${LDADD}
ORIG_LDFLAGS=${LDFLAGS}
ORIG_CXXFLAGS=${CXXFLAGS}

# Check for C compiler
AC_PROG_CXX
# Check for C++ compiler
AC_PROG_CPP

# Tells automake to create a Makefile
# See https://www.gnu.org/software/automake/manual/html_node/Requirements.html


#dnl add option for compiling with OpenMP 
AC_ARG_ENABLE([OpenMP],               
	      [AS_HELP_STRING([--disable-OpenMP],                               
			      [Disable the use of OpenMP])],               
		[],               
		[AX_OPENMP]) 

#dnl Checks for libraries.
CXXFLAGS="-Wall -std=c++11 -pedantic -funroll-loops -lz -fopenmp"
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}" 
CFLAGS="${OPENMP_CFLAGS}" 
LDFLAGS="${OPENMP_LDFLAG}"

dnl Add option --enable-debug (default to --disable-debug)
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Compile without optimization and with debug flag. Passing argument 'full' enable debugging informations.])])

AS_IF([test -z "${enable_debug}" -o "x${enable_debug}" = "xno"],
      [AC_DEFINE([NDEBUG],[],[Remove assertion checking])
       CXXFLAGS="${CXXFLAGS} -O3"],
      [CXXFLAGS="${CXXFLAGS} -g -O0 -fno-inline"
       AS_IF([test "x${enable_debug}" == "xfull"],
             [AC_DEFINE([DEBUG],[],[Output (lots of) debugging informations])]
            )]
     )

dnl =============================== Compressinon ===============================

Compressinon_LDADD=
Compressinon_LDFLAGS=
Compressinon_CXXFLAGS=
dnl Add option --without-included-Compressinon (default --with-included-Compressinon=check)
AC_ARG_WITH([included-Compressinon],
            [AS_HELP_STRING([--without-included-Compressinon],
                            [don't use the Compressinon library included here])],
            [],
            [with_included_Compressinon="check"])

dnl Add option --with-libCompressinon-prefix=PATH (default --without-libCompressinon-prefix)
AC_ARG_WITH([libCompressinon-prefix],
            [AS_HELP_STRING([--with-libCompressinon-prefix=DIR],
                            [search for libCompressinon in DIR/include/libCompressinon and DIR/lib])],
            [Compressinon_LDFLAGS="-L${with_libCompressinon_prefix}/lib"
             Compressinon_CXXFLAGS="-I${with_libCompressinon_prefix}/include/src"],
            [with_libCompressinon_prefix=""])

dnl Check if we need included Compressinon library.
PB_OK=0
AS_IF([test "x${with_included_Compressinon}" == "xcheck"],
      [AS_IF([test "x${with_libCompressinon_prefix}" != "x"],
             [with_included_Compressinon="no"],
             [with_included_Compressinon="yes"])])

OLD_LDADD=${LDADD}
OLD_LDFLAGS=${LDFLAGS}
OLD_CXXFLAGS=${CXXFLAGS}
AS_IF([test "x${with_included_Compressinon}" == "xyes"],
      [AS_IF([test "x${with_libCompressinon_prefix}" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-Compressinon and --with-libCompressinon-prefix options])],
             [PB_OK=1
              Compressinon_LDADD="\$(abs_top_builddir)/libCompressinon/src/libCompressinon.la"
              Compressinon_LDFLAGS="${Compressinon_LDFLAGS} -L\$(abs_top_builddir)/libCompressinon/src/.libs"
              Compressinon_CXXFLAGS="${Compressinon_CXXFLAGS} -I\$(abs_top_srcdir)/libCompressinon/src/"
              AC_CONFIG_SUBDIRS([libCompressinon])])],
      [AS_IF([test "x${with_libCompressinon_prefix}" != "x"],
             [
	      Compressinon_LDADD="${with_libCompressinon_prefix}/lib/libCompressinon.la"
              Compressinon_LDFLAGS="${Compressinon_LDFLAGS} -L${with_libCompressinon_prefix}/lib"
              Compressinon_CXXFLAGS="${Compressinon_CXXFLAGS} -I${with_libCompressinon_prefix}/include/libCompressinon"
	      LDADD="${LDADD} ${Compressinon_LDADD}"
	      LDFLAGS="${LDFLAGS} ${Compressinon_LDFLAGS}"
	      CXXFLAGS="${CXXFLAGS} ${Compressinon_CXXFLAGS}"
	      AC_CHECK_LIB([Compressinon],
                           [libCompressinonVersion],
                           [PB_OK=1])],
             [AC_CHECK_LIB([Compressinon],
                           [libPCompressinonVersion],
                           [PB_OK=1])])]) dnl Override the default behavior of AC_CHECK_LIB (we don't want -lCompressinon in LIBS).

LDADD=${OLD_LDADD}
LDFLAGS=${OLD_LDFLAGS}
CXXFLAGS=${OLD_CXXFLAGS}
AS_IF([test "${PB_OK}" == "1"],
      [Compressinon_LDFLAGS="${Compressinon_LDFLAGS} -lCompressinon"
       AC_DEFINE([HAVE_LIBCOMPRESSION], [1], [Define to 1 in order to Use libCompressinon])],
      [AC_MSG_FAILURE([
Unable to find a preinstalled Compressinon library.
You should run the configure script again with
'--with-included-Compressinon' option.
])])

dnl Set INCLUDED_COMPRESSION conditionnal for automake
AM_CONDITIONAL([INCLUDED_COMPRESSION], [test "x${with_included_Compressinon}" == "xyes"])

dnl Set Compressinon_{LDADD, LDFLAGS, CXXFlAGS} substituable variables
AC_SUBST([Compressinon_LDADD])
AC_SUBST([Compressinon_LDFLAGS])
AC_SUBST([Compressinon_CXXFLAGS])

dnl Forward options given to current configure script to libCompressinon configure script if needed
AS_IF([test "x${enable_shared}" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-shared"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-shared"])])
AS_IF([test "x${enable_static}" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-static"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-static"])])

dnl Forward FLAGS to libCompressinon configure scripts if needed
AS_VAR_APPEND([ac_configure_args], [" LDADD=\"${LDADD}\""])
AS_VAR_APPEND([ac_configure_args], [" LDFLAGS=\"${LDFLAGS}\""])
AS_VAR_APPEND([ac_configure_args], [" CXXFLAGS=\"${CXXFLAGS}\""])

dnl ============================ End of Compressinon ===========================


#dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([\
    libintl.h stdint.h \
    sys/ioctl.h sys/time.h \
    unistd.h zlib.h\
  ], [], [AC_MSG_ERROR([Header file '${ac_header}' not found.])])
#dnl Checks for library functions.
AC_SEARCH_LIBS([deflate], [z], [], [AC_MSG_ERROR([Library 'zlib' not found.])])

#dnl Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

#dnl Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([\
  atexit memchr memset setlocale \
  strcasecmp strdup strerror strtoul \
])

#dnl Add original existing environment variables at the end of configures variables (allow override)
LDADD="${LDADD} ${ORIG_LDADD}"
LDFLAGS="${LDFLAGS} ${ORIG_LDFLAGS}"
CXXFLAGS="${CXXFLAGS} ${ORIG_CXXFLAGS}"


AC_CONFIG_FILES([
  Makefile src/Makefile
])
AC_OUTPUT
