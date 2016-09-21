AC_INIT([flimanalysis], [0.1], [matthias.reis@biologie.hu-berlin.de])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

if test x"$debug" = x"true"; then
	CFLAGS='-O0'
else
	CFLAGS='-O2'
fi
AC_PROG_CC


# Checks for libraries.
AC_CHECK_LIB(m,sin)

AC_CHECK_HEADERS([tiffio.h])
AC_CHECK_HEADERS([tif_dir.h], [], [], [#include <tiffio.h>])

AC_CHECK_LIB([tiff], [TIFFOpen], [TIFF_LIBS="-ltiff"], AC_MSG_ERROR("Cannot build without libtiff [does your system require a libtiff-devel package?]"), -lm)
AC_SUBST(TIFF_CFLAGS)
AC_SUBST(TIFF_LIBS)

#AC_CHECK_LIB([gsl], [GSL_LIB="-lgsl", AC_MSG_ERROR("Cannot build without libgsl [does your system require a libgsl-devel package?)"), -lm)

PKG_CHECK_MODULES(GSL, gsl)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

#AC_CHECK_LIB([gslcblas], [GSL_LIB="-lgslcblas", AC_MSG_ERROR("Cannot build without libgslcblas (does your system require a libgslcblas-devel package?)"), -lm)
#PKG_CHECK_MODULES(GSLCBLAS, gslcblas)
AC_CHECK_LIB([gslcblas], [cblas_dgemm], ,[AC_MSG_ERROR(libgslcblas not found or out of sync)])
AC_SUBST(GSLCBLAS_CFLAGS)
AC_SUBST(GSLCBLAS_LIBS)


PKG_CHECK_MODULES(GLIB, glib-2.0)
AC_SUBST(GLIB_CFLAGS)
AC_SUBST(GLIB_LIBS)

CFLAGS="$CFLAGS -Wall"
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

