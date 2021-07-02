
/* model function for several ROIs */
void
multi_convolve (double *new_par, double *function, int sum_npvar, int sum_nsamples, void *data);

void convolved_exp_f_gs (double *new_par, double *function, int npvar,
                         int nsamples, void *data);

void convolved_exp_df_gs (double *new_par, double *jac, int npar,
                          int nsamples, void *data);
