#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_vector.h>



#include "fit.h"
#include "convolution.h"
#include "em.h"
#include "utils.h"
#include "levmar_mle_diff.h"
#include "fileio.h"

//#define DEBUG

void
jac_log_exp (void (*func) (double *p, double *hx, int m, int n, void *adata),
	     /* function to differentiate */
	     double *p,		/* I: current parameter estimate, mx1 */
	     double *hx,	/* I: func evaluated at p, i.e. hx=func(p), nx1 */
	     double *hxx,	/* W/O: work array for evaluating func(p+delta), nx1 */
	     double delta,	/* increment for computing the Jacobian */
	     double *jac,	/* O: array for storing approximated Jacobian, nxm */
	     int m, int n, void *adata);

void
jac_log_exp (void (*func) (double *p, double *hx, int m, int n, void *adata),
	     /* function to differentiate */
	     double *p,		/* I: current parameter estimate, mx1 */
	     double *hx,	/* I: func evaluated at p, i.e. hx=func(p), nx1 */
	     double *hxx,	/* W/O: work array for evaluating func(p+delta), nx1 */
	     double delta,	/* increment for computing the Jacobian */
	     double *jac,	/* O: array for storing approximated Jacobian, nxm */
	     int m, int n, void *adata)
{
  register int i, j;
  double tmp;
  register double d;

  for (j = 0; j < m; ++j)
    {
      (*func) (p, hx, m, n, adata);
      /* determine d=max(1E-04*|p[j]|, delta), see HZ */
      d = 1E-4 * p[j];		// force evaluation
      d = FABS (d);
      if (d < delta)
	d = delta;
      tmp = p[j];
      p[j] += d;

      (*func) (p, hxx, m, n, adata);

      printf ("hx[0]: %f, hxx[0]: %f\n", hx[0], hxx[0]);

      p[j] = tmp;		/* restore */

      d = (1.0) / d;		/* invert so that divisions can be carried out faster as multiplications */
      for (i = 0; i < n; ++i)
	{
	  jac[i * m + j] = (log (hxx[i]) - log (hx[i])) * d;
	  if (jac[i * m + j] != jac[i * m + j])
	    jac[i * m + j] = 0.0;
	}
    }
}


void
free_convolution_par (Model * conv_mdl)
{
//    UINT i;
/*
    free (conv_mdl->pvar);
    for (i = 0; i < conv_mdl->nexp; i++) free (conv_mdl->pdf[i]);
    free(conv_mdl->pdf);
*/
  free (conv_mdl->irf);
  free (conv_mdl->shifted_irf);
  free (conv_mdl->darknoise);
//    free(conv_mdl->par);
  free (conv_mdl);

  return;
}


void
free_measurement (Histogram * measure)
{
  free (measure->sample);

  free (measure);

  return;
}

/* a wrapper around dlevmar_mle_diff (Extended Levenberg-Marquardt Algorithm by Lawrence et al.
                                      using the convolved exponential model defined in convolution.c*/
uint8_t
fit_levmar_histogram (double *histogram, Model * mdl, uint8_t fitType)
{
  int status, i;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const int maxIterations = 1000;
  double *dbl_pvar = malloc (sizeof (double) * mdl->npvar);

  printf ("Beginning fit...\n");

  opts[0] = LM_INIT_MU;
  opts[1] = LM_STOP_THRESH;
  opts[2] = LM_STOP_THRESH;	//1e-200;
  opts[3] = LM_STOP_THRESH;
  opts[4] = 1e-200;

  for (i = 0; i < mdl->npvar; i++)
    {
      dbl_pvar[i] = mdl->par[mdl->pvar[i]];
    }

  status = dlevmar_mle_diff (&convolved_exp_f_gs, dbl_pvar,
			     histogram, mdl->npvar, mdl->nsamples,
			     maxIterations, opts, info, NULL, NULL, mdl,
			     fitType);


  for (i = 0; i < mdl->npvar; i++)
    {
      printf ("par[mdl->pvar[%d]]: %f, mdl->pvar[%d]: %d\n", i,
	      mdl->par[mdl->pvar[i]] = dbl_pvar[i], i, mdl->pvar[i]);
    }

  printf ("Chi square: %f\nIterations: %f\n",
	  info[1] / (mdl->tstop - mdl->tstart - mdl->npvar), info[5]);
  printf ("Exit code: %f (%s)\n", info[6], lm_error_str (info[6]));

  free (dbl_pvar);
  return info[6];
}

/* similar to fit_levmar_histogram, except that a model function for several ROIs is used*/
uint8_t
fit_levmar_FLIM_measurement (FLIM_measurement * m)
{
  int status;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const int maxIterations = 1000;

  printf ("Beginning fit...\n");

  opts[0] = LM_INIT_MU;
  opts[1] = LM_STOP_THRESH;
  opts[2] = LM_STOP_THRESH;	//1e-200;
  opts[3] = LM_STOP_THRESH;
  opts[4] = 1e-200;

#ifdef DEBUG
  plot_double (m->joint_histogram, m->nsamples);
#endif

  status = dlevmar_mle_diff (&multi_convolve, m->global_mdl->pvar_multi,
			     m->joint_histogram, m->global_mdl->npvar,
			     m->nsamples * m->nrois, maxIterations, opts,
			     info, NULL, NULL, m, m->global_mdl->fitType);

  write_levmar_result (m, info);

  return info[6];
}



uint8_t
fit_exp_gs (Histogram * measure, Model * mdl, uint8_t fitType)
{
  int status, i;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const int maxIterations = 1000;
  double *dbl_pvar = malloc (sizeof (double) * mdl->npvar);
  double e[4096], chisq, averagetau;
  FILE *fp;

  printf ("Beginning fit...\n");

  opts[0] = LM_INIT_MU;
  opts[1] = LM_STOP_THRESH;
  opts[2] = LM_STOP_THRESH;	//FLT_MIN;
  opts[3] = LM_STOP_THRESH;
  opts[4] = LM_STOP_THRESH;

  for (i = 0; i < mdl->npvar; i++)
    {
      dbl_pvar[i] = mdl->par[mdl->pvar[i]];
    }
  mdl->shifted_irf = malloc (sizeof (double) * mdl->nsamples);
  memcpy (mdl->shifted_irf, mdl->irf, mdl->nsamples * sizeof (double));
  mdl->shifted_irf_offset = 0.0;

  if (mdl->npvar != 0)
    {
      status = dlevmar_mle_diff (&convolved_exp_f_gs, dbl_pvar,
				 measure->sample, mdl->npvar,
				 measure->nsamples, maxIterations, opts, info,
				 NULL, NULL, mdl, fitType);
/*		status = dlevmar_dif( &convolved_exp_f_gs, dbl_pvar,
			measure->sample[0], mdl->npvar, measure->nsamples,
			maxIterations,opts,info,NULL,NULL,mdl);*/

      for (i = 0; i < mdl->npvar; i++)
	{
	  printf ("par[mdl->pvar[%d]]: %f, mdl->pvar[%d]: %d\n", i,
		  mdl->par[mdl->pvar[i]] = dbl_pvar[i], i, mdl->pvar[i]);
	}
      averagetau = 0.0;
      for (i = IDX_TAU0; i < mdl->nexp * 2 + IDX_TAU0; i += 2)
	{
	  averagetau += mdl->par[i] * mdl->par[i + 1];
	  printf ("Amp[%d] %f\tLifetime[%d] %f\n", i, mdl->par[i + 1], i,
		  mdl->par[i] * measure->tcal);
	}
      printf ("Average Lifetime: %f (tcal %f)\n", averagetau * measure->tcal,
	      measure->tcal);

      printf ("Chi square: %f\nIterations: %f\n",
	      info[1] / (mdl->tstop - mdl->tstart - mdl->npvar), info[5]);
      fp = fopen ("averagetau.csv", "at");
      fprintf (fp, "%f\n", averagetau);
      fclose (fp);
    }
  mdl->complete_pdf = malloc (sizeof (double) * mdl->nsamples);
  convolved_exp_f_gs (dbl_pvar, mdl->complete_pdf, mdl->npvar,
		      mdl->nsamples, mdl);
/*	plot_double (mdl->complete_pdf, mdl->nsamples);
	plot_double (measure->sample[0], mdl->nsamples);*/
  chisq =
    compute_chisq_measure (e, measure->sample, mdl->complete_pdf,
			   mdl->nsamples, fitType, NULL);
  printf ("Chi square: %f, ncounts %f\n",
	  chisq / (mdl->tstop - mdl->tstart - mdl->npvar),
	  calc_ncounts (mdl->complete_pdf, mdl->nsamples));

  free (dbl_pvar);
  return info[6];
}


char *
emerr2str (uint8_t res)
{
  char *str = malloc (sizeof (char) * 255);

  switch (res)
    {
    case MAX_IT:
      str = "Maximum number of iterations reached.";
      break;
    case SMALL_DELTA:
      str = "Smallest allowed step size reached.";
      break;
    }

  return str;
}

void
fit_FLIM_measurement (FLIM_measurement * m)
{
  uint32_t i, j;
  double ***matrix, **average_matrix, roiamp[m->global_mdl->nexp];

  matrix = malloc (m->width * sizeof (double **));
  average_matrix = malloc (m->width * sizeof (double *));
  for (i = 0; i < m->width; i++)
  {
	  matrix[i] = malloc (m->width * sizeof (double *));
	  average_matrix[i] = calloc (m->height, sizeof (double));
  }
  for (i = 0; i < m->width; i++)
  {
	  for (j = 0; j < m->height; j++)
	  {
		  matrix[i][j] = calloc (m->global_mdl->nexp, sizeof (double));
	  }
  }

  fit_levmar_FLIM_measurement (m);
  for (i = 0; i < m->nrois; i++)
  {
	  double average;
	  uint8_t comp;

	  printf ("Fitting lifetimes pixel wise for ROI %d...\n", i);

	  for (comp = 0; comp < m->global_mdl->nexp; comp++)
		  roiamp[comp] = m->roi[i].mdl->par[IDX_TAU0 + 2 * comp + 1];

	  for (j = 0; j < m->roi[i].npixels; j++)
	  {
		  if (m->roi[i].pixel[j].ncounts >= m->global_mdl->threshold)
		  {
			  fit_amplitude_tac (m->roi[i].pixel[j].photons,
					  m->roi[i].pixel[j].ncounts, m->roi[i].mdl);

			  average = 0.0;
			  for (comp = 0; comp < m->global_mdl->nexp; comp++)
			  {
				  matrix[m->roi[i].pixel[j].x][m->roi[i].pixel[j].y][comp] =
						  m->roi[i].mdl->par[IDX_TAU0 + 2 * comp + 1];
				  average += m->roi[i].mdl->par[IDX_TAU0 + 2 * comp + 1]*m->tcal * m->roi[i].mdl->par[IDX_TAU0 + 2 * comp];

				  m->roi[i].mdl->par[IDX_TAU0 + 2 * comp + 1] = roiamp[comp];
			  }
			  average_matrix[m->roi[i].pixel[j].x][m->roi[i].pixel[j].y] = average;
		  }
	  }
	  printf ("done.\n");
  }

  // TODO: handle overlapping ROIs
  printf ("Writing results to TIFF files...");
  write_flim_image (matrix, m);
  write_flim_image_average_tau (average_matrix, m);
  printf ("done.\n");

  free (matrix);
  free (average_matrix);
  return;
}


uint16_t
fit_amplitude_histogram (double *histogram, double ncounts, Model * conv_mdl)
{
  uint16_t i, res, paridx[conv_mdl->nexp];
  uint8_t comp;
  double last_amp[conv_mdl->nexp - 1];
  double delta;
  const uint16_t maxit = 1000;

  for (comp = 0; comp < conv_mdl->nexp; comp++)
    {
      paridx[comp] = IDX_TAU0 + 2 * comp + 1;
      last_amp[comp] = conv_mdl->par[paridx[comp]];
      if (last_amp[comp] == 1.0)
	{
	  conv_mdl->par[paridx[comp]] = last_amp[comp] = 0.99;	//TODO: ist das eine gute Lösung?
	}
    }

  res = MAX_IT;
  for (i = 0; i < maxit; i++)
    {
      iterate_amplitude_histogram (histogram, ncounts, conv_mdl);
      delta = 0.0;
      for (comp = 0; comp < (conv_mdl->nexp - 1); comp++)
	{
	  delta += fabs (last_amp[comp] - conv_mdl->par[paridx[comp]]);
	  last_amp[comp] = conv_mdl->par[paridx[comp]];
	}
      if (delta < 0.0001)
	{
	  res = SMALL_DELTA;
	  break;
	}
    }

  return res;
}

uint16_t
fit_amplitude_tac (GSList * photons, uint64_t ncounts, Model * conv_mdl)
{
  uint16_t i, res, paridx[conv_mdl->nexp];
  uint8_t comp;
  double last_amp[conv_mdl->nexp - 1];
  double delta;
  const uint16_t maxit = 1000;


  for (comp = 0; comp < conv_mdl->nexp; comp++)
  {
	  paridx[comp] = IDX_TAU0 + 2 * comp + 1;
	  last_amp[comp] = conv_mdl->par[paridx[comp]];
	  if (last_amp[comp] == 1.0)
	  {
		  conv_mdl->par[paridx[comp]] = last_amp[comp] = 0.99;	//TODO: ist das eine gute Lösung?
	  }
  }

  res = MAX_IT;
  for (i = 0; i < maxit; i++)
  {
	  iterate_amplitude_tac (photons, ncounts, conv_mdl);
	  delta = 0.0;

	  // Normalization
	  for (comp = 0; comp < (conv_mdl->nexp - 1); comp++)
	  {
		  delta += fabs (last_amp[comp] - conv_mdl->par[paridx[comp]]);
		  last_amp[comp] = conv_mdl->par[paridx[comp]];
	  }
	  if (delta < 0.0001)
	  {
		  res = SMALL_DELTA;
		  break;
	  }
  }

  return res;
}

/* Not tested for more than one parameter */
double *
cramer_rao_boundary_gs (Model * mdl)
{
  uint64_t i, j, k;			/* i is the timechannel, j and k are the indicies of the fisher matrix */
  int s;
  double *jac, *variance, ncounts_bak;
  gsl_matrix *fisher_information = gsl_matrix_alloc (mdl->npvar, mdl->npvar);
  gsl_matrix *inv_fisher_information =
    gsl_matrix_alloc (mdl->npvar, mdl->npvar);
  gsl_matrix_view LUdecomp_view;
  gsl_permutation *perm = gsl_permutation_alloc (mdl->npvar);
  double *dbl_pvar = malloc (sizeof (double) * mdl->npvar);
  double *hxx = malloc (sizeof (double) * mdl->nsamples);
  double *hx = malloc (sizeof (double) * mdl->nsamples);

  jac = malloc (sizeof (double) * mdl->nsamples * mdl->npvar);


  for (i = 0; i < mdl->npvar; i++)
    {
      dbl_pvar[i] = mdl->par[mdl->pvar[i]];
    }

  ncounts_bak = mdl->ncounts;
  mdl->ncounts = 1.0;
  jac_log_exp (&convolved_exp_f_gs,
	       /* function to differentiate */
	       dbl_pvar,	/* I: current parameter estimate, mx1 */
	       hx,		/* I: func evaluated at p, i.e. hx=func(p), nx1 */
	       hxx,		/* W/O: work array for evaluating func(p+delta), nx1 */
	       1e-6,		/* increment for computing the Jacobian */
	       jac,		/* O: array for storing approximated Jacobian, nxm */
	       mdl->npvar, mdl->nsamples, mdl);
  LUdecomp_view = gsl_matrix_view_array (jac, mdl->npvar, mdl->nsamples);
  mdl->ncounts = ncounts_bak;

  for (j = 0; j < mdl->npvar; j++)
    {
      for (k = 0; k < mdl->npvar; k++)
	{
	  double expected = 0.0;
	  gsl_matrix_set (fisher_information, j, k, 0.0);
	  for (i = 0; i < mdl->nsamples; i++)
	    {
	      expected += hx[i] *
		gsl_matrix_get (&LUdecomp_view.matrix, j, i) *
		gsl_matrix_get (&LUdecomp_view.matrix, k, i);
//                                      expected += mdl->sample[0][0] *exp(-((double) i)/1000.0) *
//                                              (((double) i)/(1000.0*1000.0)-1.0/1000.0)
//                                              *(((double) i)/(1000.0*1000.0)-1.0/1000.0);
	    }
	  gsl_matrix_set (fisher_information, j, k, expected);
	}
    }
  gsl_linalg_LU_decomp (fisher_information, perm, &s);
  gsl_linalg_LU_invert (fisher_information, perm, inv_fisher_information);

  variance = malloc (sizeof (double) * mdl->npvar);
  for (j = 0; j < mdl->npvar; j++)
    {
      printf ("variance: %f\n", variance[j] =
	      gsl_matrix_get (inv_fisher_information, j, j));
    }

  free (jac);
  free (hxx);
  free (hx);
  free (dbl_pvar);
  gsl_matrix_free (inv_fisher_information);
  gsl_matrix_free (fisher_information);
  gsl_permutation_free (perm);
  return variance;
}

ROI *
create_rois (uint16_t nrois, uint16_t nsamples, uint8_t nexp)
{
  uint64_t i, j;
  ROI *roi = calloc (nrois, sizeof (ROI));

  for (i = 0; i < nrois; i++)
    {
      roi[i].mdl = malloc (sizeof (Model));
      roi[i].mdl->ncounts = 0;
      roi[i].ncounts = 0;
      roi[i].npixels = 0;
      roi[i].mdl->nexp = nexp;
      roi[i].mdl->par = calloc (2*nexp + IDX_TAU0 + 1, sizeof (double));
      roi[i].mdl->npar = 2*nexp + IDX_TAU0 + 1;
      roi[i].mdl->npvar = 0;
      roi[i].mdl->pvar = calloc (2*nexp + IDX_TAU0 + 1, sizeof (uint64_t));
      roi[i].mdl->pdf = malloc (sizeof (double *)*nexp);
      for (j = 0; j < nexp; j++)
	{
	  roi[i].mdl->pdf[j] = calloc (nsamples, sizeof (double));
	}
      roi[i].mdl->complete_pdf = calloc (nsamples, sizeof (double));
      roi[i].mdl->shifted_irf = calloc (nsamples, sizeof (double));
      roi[i].mdl->nsamples = nsamples;
      roi[i].mdl->darknoise = calloc (nsamples, sizeof (double));
    }

  return roi;
}
