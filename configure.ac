#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([hvs], [0.1], [vitost@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT
#AC_CONFIG_SRCDIR([src/libhvs/hevosol.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB(m, exp)
AC_CHECK_LIB(m, pow)
AC_CHECK_LIB(m, sqrt)

# Checks for header files.
AC_CHECK_HEADERS([stdio.h string.h stdlib.h math.h unistd.h float.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset memcpy getopt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
