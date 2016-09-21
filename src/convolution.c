#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fit.h"
#include "convolution.h"
#include "utils.h"
#include "levmar_mle_diff.h"	// for fitType definition



/* model function for several ROIs */
void
multi_convolve (double *new_par, double *function, int sum_npvar,
		int sum_nsamples, void *data)
{
  FLIM_measurement *m = (FLIM_measurement *) data;
  uint64_t i, iroi, ipar, ipvar, ilink, j;
  Reference *ref;
  GList *ref_list;
  double *init_function = function;

  // update each ROIs parameters (m->roi[iroi].mdl->par[ipar]) from new_par
  for (i = 0; i < m->global_mdl->npvar; i++)
    {
//              m->global_mdl->pvar_multi[i] = new_par[i];

      for (ref_list =
	   g_list_first (((GList*) g_list_nth(m->global_mdl->pvar_multi_map,i))->data);
	   ref_list; ref_list = g_list_next (ref_list))
	{

	  ref = (Reference *) ref_list->data;
	  ipar = ref->ipar;
	  ipvar = ref->ipvar;
	  iroi = ref->iroi;
	  ilink = ref->ilink;

	  m->roi[iroi].mdl->par[ipar] = new_par[i];
	  if (ilink)
	    {
	      m->global_mdl->link[ilink - 1] = new_par[i];
	    }
	}
    }

  for (i = 0; i < m->nrois; i++)
    {
      double dbl_pvar[m->roi[i].mdl->npvar];
//              printf ("roi: %d\n", i);
      for (j = 0; j < m->roi[i].mdl->npvar; j++)
	{
	  dbl_pvar[j] = m->roi[i].mdl->par[m->roi[i].mdl->pvar[j]];
	}
      convolved_exp_f_gs (dbl_pvar, function, m->roi[i].mdl->npvar,
			  m->roi[i].mdl->nsamples, (void *) m->roi[i].mdl);

//              plot_double (function, m->nsamples);
      memset (&(function[m->roi[i].mdl->tstop]), 0.0,
	      (m->nsamples - m->roi[i].mdl->tstop) * sizeof (double));
      memset (function, 0.0, m->roi[i].mdl->tstart * sizeof (double));
      function += m->roi[i].mdl->nsamples;	// m->nsamples;
    }

  function = init_function;
//      plot_double (function, m->nsamples*m->nrois);

  return;
}


void
convolved_exp_f_gs (double *new_par, double *function, int npvar,
		    int nsamples, void *data)
{
  Model *conv_mdl = (Model *) data;

  double *irf = conv_mdl->irf;
  double *par = conv_mdl->par;
  uint16_t *pvar = conv_mdl->pvar;
  double **pdf = conv_mdl->pdf;
  double *darknoise = conv_mdl->darknoise;

  int i, j;

  double tau, amp, exp_factor[conv_mdl->nexp], trapez_factor[conv_mdl->nexp],
    partial_exp[conv_mdl->nexp];
  double sum_amp;

  int numExponentials = conv_mdl->nexp;

//      if ( (nsamples != mdl->n ) || (npar != mdl->m) ) return;

//    printf ("Calculating function...\n");

  /* Write the variable parameters into the array par */
  /* TODO: überflüssig, da bereits in multi_conv erledigt */
  for (i = 0; i < npvar; i++)
    par[pvar[i]] = new_par[i];
  /* Calculate last amplitude */
  if (conv_mdl->fitType == LM_CHISQ_MLE)
    {
      sum_amp = 0.0;
      for (i = 0; i < numExponentials - 1; i++)
	{			//IDX_TAU0 + 1; i < conv_mdl->npar - 2; i += 2) {
	  sum_amp += par[2 * i + IDX_TAU0 + 1];	//par[i];
	}
      par[2 * i + IDX_TAU0 + 1] = 1.0 - sum_amp;
    }
//    exp_factor = malloc(sizeof(double) * numExponentials);
//    trapez_factor = malloc(sizeof(double) * numExponentials);
//    partial_exp = malloc(sizeof(double) * numExponentials);

  for (j = 0; j < numExponentials; j++)
    {
      amp = par[2 * j + IDX_TAU0 + 1];
      tau = par[2 * j + IDX_TAU0];

//      printf("amp[%d]: %f\n", j, amp);
//      printf("tau[%d]: %f\n", j, tau*10.9);

/*	if ((tau != tau) || (tau <= 0.0)) {
		memset(function, 0, nsamples * sizeof(double));
		for (j = 0; j < numExponentials; j++) {
			memset(conv_mdl->pdf[j], 0, nsamples * sizeof(double));
		}
		conv_mdl->complete_pdf = function;
	}*/

      exp_factor[j] = exp (-1.0 / tau);
      trapez_factor[j] = 0.5;	// * amp;
      partial_exp[j] = 0.0;
    }

  if (par[IDX_IRFSHIFT] != 0.0)
    {
      memcpy (conv_mdl->shifted_irf, conv_mdl->irf,
	      conv_mdl->nsamples * sizeof (double));
      shift_histogram (conv_mdl->shifted_irf, nsamples, par[IDX_IRFSHIFT]);	//conv_mdl->shifted_irf_offset);
      conv_mdl->shifted_irf_offset = par[IDX_IRFSHIFT];
      norm_histogram (conv_mdl->shifted_irf, nsamples, 1.0);
      irf = conv_mdl->shifted_irf;
    }
  memset (function, 0.0, nsamples * sizeof (double));

  for (j = 0; j < numExponentials; j++)
    {
//      memset(pdf[j], 0.0, nsamples * sizeof(double));
      partial_exp[j] = pdf[j][0] = trapez_factor[j] * irf[0];
//      function[0] += pdf[j][0];
      for (i = 1; i < nsamples; i++)
	{
	  partial_exp[j] =
	    (partial_exp[j] +
	     trapez_factor[j] * irf[i - 1]) * exp_factor[j] +
	    trapez_factor[j] * irf[i];
	  pdf[j][i] = partial_exp[j];
//          function[i] += partial_exp[j];
	}

      norm_histogram (pdf[j], nsamples, 1.0);
      for (i = 0; i < nsamples; i++)
	{
	  function[i] += par[2 * j + IDX_TAU0 + 1] * pdf[j][i];
	}
    }
/*
    for (j = 0; j < numExponentials; j++) {
	printf ("ugu%f\n", conv_mdl->norm_factor = norm_histogram (pdf[j], nsamples,
		conv_mdl->ncounts * par[2 * j + IDX_TAU0 + 1]));
    }
*/

  conv_mdl->complete_pdf = function;
  if (conv_mdl->fitType == LM_CHISQ_MLE)
    conv_mdl->norm_factor =
      norm_histogram (function, nsamples, conv_mdl->ncounts);

// TODO: normierung und background ????
  for (i = 0; i < nsamples; i++)
    {
      function[i] += par[IDX_BG] * darknoise[i];
    }

  return;
}

void
convolved_exp_df_gs (double *new_par, double *jac, int npar, int nsamples,
		     void *data)
{
  Model *conv_mdl = (Model *) data;
  double *irf = conv_mdl->irf;
  double *par = conv_mdl->par;
  uint16_t *pvar = conv_mdl->pvar;
  double **pdf = conv_mdl->pdf;
  double *jac_row;		/* rows are for parameters, columns for time channels */

  int i, j, k;

  double tau, amp, *exp_factor, *trapez_factor, one_div_tau2;
  double sum_amp;

  int numExponentials = conv_mdl->nexp;
  int exp_idx;

  printf ("calculating derivative...\n");

//      if ( (nsamples != mdl->n ) || (npar != mdl->m) ) return;

  /* Write the parameters which are going to be fitted into the array a */
  for (i = 0; i < npar; i++)
    par[pvar[i]] = new_par[i];

  /* Calculate last amplitude */
  sum_amp = 0.0;
  for (i = 0; i < numExponentials - 1; i++)
    {
      sum_amp += par[2 * i + IDX_TAU0 + 1];
    }
  par[2 * numExponentials + IDX_TAU0 + 1] = 1.0 - sum_amp;

  exp_factor = malloc (sizeof (double) * numExponentials);
  trapez_factor = malloc (sizeof (double) * numExponentials);

  for (j = 0; j < numExponentials; j++)
    {
      amp = par[2 * j + IDX_TAU0 + 1];
      tau = par[2 * j + IDX_TAU0];

      exp_factor[j] = exp (-1.0 / tau);
      trapez_factor[j] = 0.5 * amp;
    }

  for (k = 0; k < npar; k++)
    {
      // TODO: indizies!
      /* Entries for background noise (Parameter 0) */
      switch (pvar[k])
	{
	case IDX_BG:
	  jac_row = jac;
	  for (i = 0; i < nsamples; i++, jac_row += npar)
	    jac_row[k] = conv_mdl->darknoise[i];
	  continue;
	  break;
	case IDX_IRFSHIFT:
	  break;		// TODO: muss ganz zum schluss gemacht werden
	default:
	  break;
	}

      exp_idx = (pvar[k] - IDX_TAU0) / 2;
      switch ((pvar[k] - IDX_TAU0) % 2)
	{
	  double test;
	case 0:		/* Lifetime */
	  one_div_tau2 = 1.0 / (par[pvar[k]] * par[pvar[k]]);

	  printf ("lifetime[%d]:%f\n", exp_idx, par[pvar[k]]);

	  jac_row = jac;
	  jac_row[k] = 0.0;
	  jac_row += npar;
	  test = 0.0;
	  for (i = 1; i < nsamples; i++, jac_row += npar)
	    {
	      test += pdf[exp_idx][i - 1];
	      //printf ("(jac_row-npar)[%d]: %f\n", k, (jac_row-npar)[k]);
	      jac_row[k] =
		((jac_row - npar)[k] +
		 one_div_tau2 *
		 (pdf[exp_idx][i - 1] +
		  trapez_factor[exp_idx] * irf[i - 1])) * exp_factor[exp_idx];
	    }
/*	    jac_row = jac + npar;

	    for (i = 1; i < nsamples; i++, jac_row += npar) {
		test += pdf[exp_idx][i];
		jac_row[k] -=  0.0009278 * pdf[exp_idx][i-1]; // 0.008 * (jac_row-npar)[k]
	    }*/
	  break;
	case 1:
	  printf ("amplitude[%d]:%f\n", exp_idx, par[pvar[k]]);
	  jac_row = jac;

	  for (i = 0; i < nsamples; i++, jac_row += npar)
	    {
	      jac_row[k] = pdf[exp_idx][i] / par[pvar[k]];
	    }
	  break;
	default:
	  break;
	}
    }

  {
    double jac_tau[4096];

    jac_row = jac;
    for (i = 0; i < nsamples; i++, jac_row += npar)
      {
	jac_tau[i] = jac_row[0];
      }

    plot_double (jac_tau, nsamples);
  }


  free (exp_factor);
  free (trapez_factor);

  return;
}
