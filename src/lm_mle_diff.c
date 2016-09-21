/////////////////////////////////////////////////////////////////////////////////
// 
//  Levenberg - Marquardt non-linear minimization algorithm
//  Modified and simplified by Ted Laurence to use for MLE of Poisson-distributed data; Used only for 
//              double precision, without constraints
//  Copyright (C) 2004  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  Modifications Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory. Written by Ted Laurence (laurence2@llnl.gov)
//  LLNL-CODE-424602 All rights reserved.
//  This file is part of dlevmar_mle_der
//
//  Please also read Our Notice and GNU General Public License.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>

#include "fit.h"
#include "levmar_mle_diff.h"
#include "lm_mle_compiler.h"
  
#define EPSILON       1E-12
#define ONE_THIRD     0.3333333334 /* 1.0/3.0 */
#define LM_REAL_MAX FLT_MAX
#define LM_REAL_MIN -FLT_MAX
    int
write_file (double *p, int n) 
{
  FILE * file;
  int i = 0;
   file = fopen ("file.txt", "w");	/* apend file (add text to
					   a file or create a file if it does not exist. */
  for (i = 0; i < n; i++)
    fprintf (file, "%f\n", p[i]);	/*writes */
  fclose (file);		/*done! */
  return 0;
}

 double
compute_chisq_measure (double *e, double *x, double *hx, int n, int fitType,
		       void *adata) 
{
  int i;
  double chisq = 0.0;
  
switch (fitType)
    
    {
    case LM_CHISQ_MLE:
      for (i = 0; i < n; i++)
	
	{
	  if (hx[i] > 0)
	    
	    {
	      if (x[i] == 0)
		chisq += 2 * hx[i];
	      
	      else
		chisq += 2 * (hx[i] - x[i] - x[i] * log (hx[i] / x[i]));
	      e[i] = x[i] / hx[i] - 1.0;
	    }
	  
	  else if (hx[i] != hx[i])
	    {
	      chisq += e[i] = 1000.0;
	    }
	  
	  else
	    e[i] = 0.0;
	}
      break;
    case LM_CHISQ_NEYMAN:
      for (i = 0; i < n; i++)
	
	{
	  if (hx[i] == hx[i])
	    {
	      if (x[i] == 0)
		
		{
		  chisq += (x[i] - hx[i]) * (x[i] - hx[i]);
		  e[i] = x[i] - hx[i];
		}
	      
	      else
		
		{
		  chisq += (x[i] - hx[i]) * (x[i] - hx[i]) / x[i];
		  e[i] = (x[i] - hx[i]) / x[i];
		}
	    }
	  else
	    {
	      chisq += e[i] = 1000.0;
	    }
	}
      break;
    case LM_CHISQ_EQUAL_WT:
      for (i = 0; i < n; i++)
	
	{
	  e[i] = x[i] - hx[i];
	  chisq += e[i] * e[i];
	}
      break;
    }
  return chisq;
}

  void
dlevmar_mle_fdif_forw_jac_approx (void (*func)
				  (double *p, double *hx, int m, int n,
				   void *adata), 
				  /* function to differentiate */ 
				  double *p, /* I: current parameter estimate, mx1 */ 
				  double *hx, /* I: func evaluated at p, i.e. hx=func(p), nx1 */ 
				  double *hxx, /* W/O: work array for evaluating func(p+delta), nx1 */ 
				  double delta, /* increment for computing the Jacobian */ 
				  double *jac, /* O: array for storing approximated Jacobian, nxm */ 
				  int m, int n, void *adata) 
{
  register int i, j;
  double tmp;
  register double d;
  
//double dif[4096];
    for (j = 0; j < m; ++j)
    {
      
	/* determine d=max(1E-04*|p[j]|, delta), see HZ */ 
	d = 1E-4 * p[j];	// force evaluation
      d = FABS (d);
      if (d < delta)
	d = delta;
      tmp = p[j];
      p[j] += d;
      (*func) (p, hxx, m, n, adata);
       p[j] = tmp;		/* restore */
       d = (1.0) / d;		/* invert so that divisions can be carried out faster as multiplications */
      for (i = 0; i < n; ++i)
	{
	  jac[i * m + j] = (hxx[i] - hx[i]) * d;
	  
//      dif[i] = hxx[i]-hx[i];
	}
      
//    plot_double (jac, n*m);
    }
}

 
/* 
 * This function seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function requires an analytic Jacobian. 
 *
 * Returns the number of iterations (>=0) if successful, LM_ERROR if failed
 *
 * For more details, see K. Madsen, H.B. Nielsen and O. Tingleff's lecture notes on 
 * non-linear least squares at http://www.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 */ 
int
  dlevmar_mle_diff (void (*func)
		    (double *p, double *hx, int m, int n, void *adata),
		    /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
		    double *p, /* I/O: initial parameter estimates. On output has the estimated solution */ 
		    double *x, /* I: measurement vector. NULL implies a zero vector */ 
		    int m, /* I: parameter vector dimension (i.e. #unknowns) */ 
		    int n, /* I: measurement vector dimension */ 
		    int itmax, /* I: maximum number of iterations */ 
		    double opts[5],	/* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3,\epsilon4]. Respectively the scale factor for initial \mu,
					 * stopping thresholds for ||J^T e||_inf, ||Dp||_2, chisq, and delta_chisq. Set to NULL for defaults to be used
					 */ 
		    double info[LM_INFO_SZ], 
		    /* O: information regarding the minimization. Set to NULL if don't care
		     * info[0]= ||e||_2 at initial p.
		     * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
		     * info[5]= # iterations,
		     * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
		     *                                 2 - stopped by small Dp
		     *                                 3 - stopped by itmax
		     *                                 4 - singular matrix. Restart from current p with increased mu 
		     *                                 5 - no further error reduction is possible. Restart with increased mu
		     *                                 6 - stopped by small ||e||_2
		     *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
		     *                                                                   8 - stopped by small change in chisq on successful iteration (dF<eps4)
		     * info[7]= # function evaluations
		     * info[8]= # Jacobian evaluations
		     * info[9]= # linear systems solved, i.e. # attempts for reducing error
		     */ 
		    double *work, /* working memory at least LM_DER_WORKSZ() reals large, allocated if NULL */ 
		    double *covar, /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */ 
		    void *adata,	/* pointer to possibly additional data, passed uninterpreted to func & jacf.
					 * Set to NULL if not needed
					 */ 
		    int fitType)	/* Which type of fit to perform
					 *  LM_CHISQ_MLE is MLE for Poisson distribution
					 *  LM_CHISQ_NEYMAN is least squares where sigma^2 is set to max(x,1)
					 *  LM_CHISQ_EQUAL_WT is least squares where sigma^2 is set to 1 */  
{
  register int i, j, k;	//, l;
  int worksz, freework = 0, issolved;
  
/* temp work arrays */ 
  double *e, /* nx1 */ 
   *e_test, /* nx1 */ 
   *hx, /* \hat{x}_i, nx1 */ 
   *jacTe, /* J^T e_i mx1 */ 
   *jac, /* nxm */ 
   *jacTjac, /* mxm */ 
   *LUdecomp, /* mxm */ 
   *Dp, /* mx1 */ 
   *diag_jacTjac, /* diagonal of J^T J, mx1 */ 
   *pDp, /* p + Dp, mx1 */ 
   *wrk;			/* \hat{xx}_i, nx1, used for forward differences */
   register double mu, /* damping constant */ 
    tmp;			/* mainly used in matrix & vector multiplications */
  double p_chisq, jacTe_inf, pDp_chisq;	/* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
  double p_L2, Dp_L2 = LM_REAL_MAX, dF, dL;
  double tau, eps1, eps2, eps2_sq, eps3, eps4;
  double init_p_chisq;
  const double min_hx = 1.0 / sqrt (LM_REAL_MAX);
  int nu = 2, nu2, stop = 0, nfev, njev = 0, nlss = 0;	//,K=(m>=10)? m: 10, updjac, updp=1, newjac;
  const int nm = n * m;
   gsl_matrix_view LUdecomp_view;
  gsl_vector_view jacTe_view, Dp_view;
  gsl_permutation * perm = gsl_permutation_alloc (m);
  int signum;
   mu = jacTe_inf = 0.0;	/* -Wall */
   if (n < m)
    {
      fprintf (stderr,
		"(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n",
		n, m);
      return LM_ERROR;
    }
   
    //if(!jacf){
    //  fprintf(stderr, "No function specified for computing the Jacobian");
    //  return LM_ERROR;
    //}
    if (opts)
    {
      tau = opts[0];
      eps1 = opts[1];
      eps2 = opts[2];
      eps2_sq = opts[2] * opts[2];
      eps3 = opts[3];
      eps4 = opts[4];
    }
  
  else
    {				// use default values
      tau = LM_INIT_MU;
      eps1 = LM_STOP_THRESH;
      eps2 = LM_STOP_THRESH;
      eps2_sq = LM_STOP_THRESH * LM_STOP_THRESH;
      eps3 = LM_STOP_THRESH;
      eps4 = LM_STOP_THRESH;
    }
   if (!work)
    {
      worksz = LM_MLE_WORKSZ (m, n);	//2*n+4*m + n*m + m*m; 
      work = (double *) malloc (worksz * sizeof (double));	/* allocate a big chunk in one step */
      if (!work)
	{
	  fprintf (stderr, "(): memory allocation request failed\n");
	  return LM_ERROR;
	}
      freework = 1;
    }
   
    /* set up work arrays */ 
    e = work;
  e_test = e + n;
  hx = e_test + n;
  wrk = hx + n;
  jacTe = wrk + n;
  jac = jacTe + m;
  jacTjac = jac + nm;
  LUdecomp = jacTjac + m * m;
  Dp = LUdecomp + m * m;
  diag_jacTjac = Dp + m;
  pDp = diag_jacTjac + m;
   LUdecomp_view = gsl_matrix_view_array (LUdecomp, m, m);
  jacTe_view = gsl_vector_view_array (jacTe, m);
  Dp_view = gsl_vector_view_array (Dp, m);
   
    /* compute e=x - f(p) and its chisq measure */ 
    (*func) (p, hx, m, n, adata);
  nfev = 1;
  
    /* ### e=x-hx, p_eL2=||e|| */ 
    p_chisq = compute_chisq_measure (e, x, hx, n, fitType, adata);
   init_p_chisq = p_chisq;
  if (!LM_FINITE (p_chisq))
    stop = 7;
   for (k = 0; k < itmax && !stop; ++k)
    {
      register int l, im;
      register double alpha, wt, *jaclm;
       
	/* Note that p and e have been updated at a previous iteration */ 
	if (p_chisq <= eps3)
	{			/* error is small */
	  stop = 6;
	  break;
	}
       
	/* Compute the Jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
	 * Since J^T J is symmetric, its computation can be sped up by computing
	 * only its upper triangular part and copying it to the lower part
	 */ 
	
	//(*jacf)(p, jac, m, n, adata); ++njev;
	//if((updp && nu>16) || updjac==K){ /* compute difference approximation to J */
	 dlevmar_mle_fdif_forw_jac_approx (func, p, hx, wrk, eps4, jac, m,
					     n, adata);
      
	//write_file(jac,n);
//      plot_double (jac, n);
	nfev += m;
      
	//nu=2;
	//}
	
	/* J^T J, J^T e */ 
	for (i = m * m; i-- > 0;)
	jacTjac[i] = 0.0;
      for (i = m; i-- > 0;)
	jacTe[i] = 0.0;
       for (l = n; l-- > 0;)
	{
	  jaclm = jac + l * m;
	  switch (fitType)
	    
	    {
	    case LM_CHISQ_MLE:
	      if (hx[l] > min_hx)
		wt = x[l] / (hx[l] * hx[l]);
	      
	      else
		wt = 0.0;
	      break;
	    case LM_CHISQ_NEYMAN:
	      if (x[l] == 0)
		wt = 1.0;
	      
	      else
		wt = 1.0 / x[l];
	      break;
	    case LM_CHISQ_EQUAL_WT:
	      wt = 1.0;
	      break;
	    }
	  for (i = m; i-- > 0;)
	    {
	      im = i * m;
	      alpha = jaclm[i];	//jac[l*m+i];
	      for (j = i + 1; j-- > 0;)	/* j<=i computes lower triangular part only */
		jacTjac[im + j] += jaclm[j] * alpha * wt;	//jac[l*m+j]
	      
		/* J^T e */ 
		jacTe[i] += alpha * e[l];
	    }
	}
       for (i = m; i-- > 0;)	/* copy to upper part */
	for (j = i + 1; j < m; ++j)
	  jacTjac[i * m + j] = jacTjac[j * m + i];
        
	/* Compute ||J^T e||_inf and ||p||^2 */ 
	for (i = 0, p_L2 = jacTe_inf = 0.0; i < m; ++i)
	{
	  if (jacTe_inf < (tmp = FABS (jacTe[i])))
	    jacTe_inf = tmp;
	   diag_jacTjac[i] = jacTjac[i * m + i];	/* save diagonal entries so that augmentation can be later canceled */
	  p_L2 += p[i] * p[i];
	}
       
	/* check for convergence */ 
	if ((jacTe_inf <= eps1))
	{
	  Dp_L2 = 0.0;		/* no increment for p in this case */
	  stop = 1;
	  break;
	}
       
	/* compute initial damping factor */ 
	if (k == 0)
	{
	  for (i = 0, tmp = LM_REAL_MIN; i < m; ++i)
	    if (diag_jacTjac[i] > tmp)
	      tmp = diag_jacTjac[i];	/* find max diagonal element */
	  mu = tau * tmp;
	}
       
	/* determine increment using adaptive damping */ 
	while (1)
	{
	  
	    /* augment normal equations */ 
	    for (i = 0; i < m; ++i)
	    jacTjac[i * m + i] += mu;
	   
	    /* solve augmented equations */ 
	    /* use the LU included with GSL */ 
	    for (i = 0; i < m * m; i++)
	    LUdecomp[i] = jacTjac[i];
	  gsl_linalg_LU_decomp (&LUdecomp_view.matrix, perm, &signum);
	  issolved =
	    !gsl_linalg_LU_solve (&LUdecomp_view.matrix, perm,
				  &jacTe_view.vector, &Dp_view.vector);
	  if (issolved)
	    {
	      
		/* compute p's new estimate and ||Dp||^2 */ 
		for (i = 0, Dp_L2 = 0.0; i < m; ++i)
		{
		  pDp[i] = p[i] + (tmp = Dp[i]);
		  Dp_L2 += tmp * tmp;
		}
	      
		//Dp_L2=sqrt(Dp_L2);
		if (Dp_L2 <= eps2_sq * p_L2)
		{		/* relative change in p is small, stop */
		  stop = 2;
		  break;
		}
	       
//       if(Dp_L2>=(p_L2+eps2)/EPSILON*EPSILON){ /* almost singular */
//         stop=4;
//         break;
//      }
		(*func) (pDp, hx, m, n, adata);
	      ++nfev;		/* evaluate function at p + Dp */
	      
		/* compute ||e(pDp)||_2 */ 
		/* ### hx=x-hx, pDp_chisq=||hx|| */ 
		pDp_chisq =
		compute_chisq_measure (e_test, x, hx, n, fitType, adata);
	       if (!LM_FINITE (pDp_chisq))
		{		/* chisq is not finite, most probably due to a user error.
				 * This check makes sure that the inner loop does not run indefinitely.
				 * Thanks to Steve Danauskas for reporting such cases
				 */
		  stop = 7;
		  break;
		}
	       for (i = 0, dL = 0.0; i < m; ++i)
		dL += Dp[i] * (mu * Dp[i] + jacTe[i]);
	       dF = p_chisq - pDp_chisq;
	       if (dL > 0.0 && dF > 0.0)
		{		/* reduction in error, increment is accepted */
		  tmp = (2.0 * dF / dL - 1.0);
		  tmp = 1.0 - tmp * tmp * tmp;
		  mu = mu * ((tmp >= ONE_THIRD) ? tmp : ONE_THIRD);
		  nu = 2;
		   for (i = 0; i < m; ++i)	/* update p's estimate */
		    p[i] = pDp[i];
		   for (i = 0; i < n; ++i)	/* update e and ||e||_2 */
		    e[i] = e_test[i];
		  p_chisq = pDp_chisq;
		  if (dF < eps4)
		    stop = 8;
		  break;
		}
	    }
	   
	    /* if this point is reached, either the linear system could not be solved or
	     * the error did not reduce; in any case, the increment must be rejected
	     */ 
	     mu *= nu;
	  nu2 = nu << 1;	// 2*nu;
	  if (nu2 <= nu)
	    {			/* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
	      stop = 5;
	      break;
	    }
	  nu = nu2;
	   for (i = 0; i < m; ++i)	/* restore diagonal J^T J entries */
	    jacTjac[i * m + i] = diag_jacTjac[i];
	}			/* inner loop */
    }
   if (k >= itmax)
    stop = 3;
   for (i = 0; i < m; ++i)	/* restore diagonal J^T J entries */
    jacTjac[i * m + i] = diag_jacTjac[i];
   if (info)
    {
      info[0] = init_p_chisq;
      info[1] = p_chisq;
      info[2] = jacTe_inf;
      info[3] = Dp_L2;
      for (i = 0, tmp = LM_REAL_MIN; i < m; ++i)
	if (tmp < jacTjac[i * m + i])
	  tmp = jacTjac[i * m + i];
      info[4] = mu / tmp;
      info[5] = (double) k;
      info[6] = (double) stop;
      info[7] = (double) nfev;
      info[8] = (double) njev;
      info[9] = (double) nlss;
    }
   if (freework)
    free (work);
   gsl_permutation_free (perm);
   return (stop != 4 && stop != 7) ? k : LM_ERROR;
}

  char *
lm_error_str (int error)
{
  char *error_str;
   error_str = malloc (sizeof (char) * 255);
   switch (error)
    {
    case 1:
      error_str = "stopped by small gradient J^T e";
      break;
    case 2:
      error_str = "stopped by small Dp";
      break;
    case 3:
      error_str = "stopped by itmax";
      break;
    case 4:
      error_str =
	"singular matrix. Restart from current p with increased mu";
      break;
    case 5:
      error_str =
	"no further error reduction is possible. Restart with increased mu";
      break;
    case 6:
      error_str = "stopped by small ||e||_2";
      break;
    case 7:
      error_str =
	"stopped by invalid (i.e. NaN or Inf) \"func\" values. This is a user error";
      break;
    case 8:
      error_str =
	"stopped by small change in chisq on successful iteration (dF<eps4)";
      break;
    default:
      error_str = "Internal error in lm_error_str ()";
      break;
    }
   return error_str;
}

char *
fitType_str (int fitType)
{
  char *lmfitType = malloc (sizeof (char) * 100);
   switch (fitType)
    {
    case LM_CHISQ_MLE:
      lmfitType = "Maximum-Likelihood for Poisson distribution";
      break;
    case LM_CHISQ_NEYMAN:
      lmfitType = "Neyman-Chisquare (sigma^2 is set to max(x,1))";
      break;
    case LM_CHISQ_EQUAL_WT:
      lmfitType = "Chisquare (sigma^2 is set to 1)";
      break;
    default:
      break;
    }
   return lmfitType;
}


