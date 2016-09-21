void write2file (double *array, uint64_t size);
void read_from_file (double *array, uint64_t size);
void write2expfile (double *amp, double *lifetime, double *meanamp,
		    uint64_t meanampsize, double *ncounts, uint64_t ndifcounts,
		    uint64_t nmodels);
void write2expampfixfile (double *variance, double variance_boundary,
			  double *ncounts, uint64_t ndifcounts, double *meanamp);

void write_flim_image (double ***amp, FLIM_measurement * m);
void write_flim_image_average_tau (double **matrix, FLIM_measurement * m);
void write_levmar_result (FLIM_measurement * m, double *info);

void read_histogram_from_file (double *histogram, uint64_t total_size,
			       char *fname, uint64_t start, uint64_t end);
void read_globals_file (Histogram * measure, Model * conv_mdl, char *fname);
FLIM_measurement *read_flim_measurement (char *ini_file);
int
calculate_roi_pixels_from_tiff_file (char *filename, FLIM_measurement *m);
