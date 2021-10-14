# Must init the autoconf setup
# The first parameter is project name
# second is version number
# third is bug report address
AC_INIT([nihm], [1.0])

# Safety checks in case user overwritten --srcdir
AC_CONFIG_SRCDIR([src/nihm.cpp])

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

${CXXFLAGS="-g -O3 -Wall -std=c++11 -pedantic -lz -fopenmp -flto -march=native -mtune=native"}
${LDFLAGS="-lz -fopenmp"}


# Check for C compiler
AC_PROG_CXX
# Check for C++ compiler
AC_PROG_CPP

# Tells automake to create a Makefile
# See https://www.gnu.org/software/automake/manual/html_node/Requirements.html
AC_CONFIG_FILES([
		 Makefile src/Makefile
		 ])

# Generate the output
AC_OUTPUT

#dnl add option for compiling with OpenMP 
AC_ARG_ENABLE([OpenMP],               
	      [AS_HELP_STRING([--disable-OpenMP],                               
			      [Disable the use of OpenMP])],               
		[],               
		[AX_OPENMP]) 
dnl Checks for libraries. 
CXXFLAGS="-Wall -std=c++11 -pedantic -funroll-loops" 
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}" 
CFLAGS="${OPENMP_CFLAGS}" 
LDFLAGS="${OPENMP_LDFLAG}"
