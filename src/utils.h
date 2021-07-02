
#include <stdint.h>

#include "gnuplot_i.h"

double norm_histogram (double *histogram, uint64_t ntchannels, double ncounts);
double norm_histogram_peak (double *histogram, uint64_t ntchannels, uint64_t ncounts);
int shift_histogram (double *a, int N, double s);

double calc_ncounts (double *histogram, uint64_t size);

void plot_double (double *sample, uint64_t ntchannels);
void plot_double_tex (double *sample, char *xtitle, char *ytitle,
		      uint64_t ntchannels);
void plot_double_xy_tex (double *x, char *xtitle, double *y, char *ytitle,
			 uint64_t ntchannels);
void plot_double_xy_yerror_tex (double *x, char *xtitle, double *y,
				char *ytitle, double *yerror,
				uint64_t ntchannels);
void setup_gnuplot_tex (gnuplot_ctrl * ploth, char *fname, FILE * gnuFile);
