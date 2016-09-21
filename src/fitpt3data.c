#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fit.h"
#include "levmar_mle_diff.h"
#include "fileio.h"
#include "utils.h"
#include "readpt3.h"
#include "ini.h"


int main(int argc, char **argv) {
	FLIM_measurement *m;

	if ((m = read_flim_measurement (argv[1])) == NULL) {
		printf ("Error while reading FLIM measurement. Exiting...\n");
		return 1;
	}
/*
	plot_double (m->joint_histogram, m->nsamples * m->nrois);
	plot_double (m->joint_histogram_irf, m->nsamples * m->nirf_rois);
*/
	fit_FLIM_measurement (m);

	return 0;
}

