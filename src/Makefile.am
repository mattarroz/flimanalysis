#AM_CPPFLAGS = $(TIFF_CFLAGS)
#pt3tiff_LDADD = $(TIFF_LIBS)
if DEBUG
AM_CFLAGS = $(GSL_CFLAGS) $(GLIB_CFLAGS) $(TIFF_CFLAGS) -ggdb -DDEBUG
else
AM_CFLAGS = $(GSL_CFLAGS) $(GLIB_CFLAGS) $(TIFF_CFLAGS)
endif

flimanalysis_LDADD = $(GSL_LIBS) $(GLIB_LIBS) $(TIFF_LIBS)

bin_PROGRAMS = flimanalysis
flimanalysis_SOURCES = pt3io.c utils.c lm_mle_diff.c fit.c convolution.c gnuplot_i.c em.c fileio.c readpt3.c read_qadata.c ini.c ferrorf.c chgext.c fitpt3data.c  fitpt3data.h utils.h  fit.h convolution.h gnuplot_i.h em.h fileio.h readpt3.h read_qadata.h ini.h errors.h pt3io.h
