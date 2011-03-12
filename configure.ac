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

AC_OPENMP

ac_openmp30=no
if test -n "${OPENMP_CFLAGS}"; then
    AC_MSG_CHECKING([for support of OpenMP 3.0])
    ac_save_CFLAGS=$CFLAGS
    CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
    AC_LINK_IFELSE( [AC_LANG_PROGRAM([
#include <omp.h>
    ],
    [[
int main() {
    int i,j,k[100];
#pragma omp parallel for collapse(2) private(i,j) shared(k)
    for(i=0;i<10;i++)
	for(j=0;j<10;j++) {
		k[i*10+j]=1;
	}
    return 0;
}
    ]])], [ac_openmp30=yes],[ac_openmp30=no] )
    CFLAGS=$ac_save_CFLAGS
    AC_MSG_RESULT([$ac_openmp30])
fi

if test "${ac_openmp30}" = "no"; then
    OPENMP_CFLAGS=""
    AC_MSG_WARN([OpenMP 3.0 is not supported])
fi

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
