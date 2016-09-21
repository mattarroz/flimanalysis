/* 
////////////////////////////////////////////////////////////////////////////////////
// 
//  Prototypes and definitions for the Levenberg - Marquardt minimization algorithm
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
////////////////////////////////////////////////////////////////////////////////////
*/
 
#include <float.h>
#ifndef _LEVMAR_H_
#define _LEVMAR_H_
   
#ifdef __cplusplus
extern "C"
{
  
#endif				/*  */
  
#define FABS(x) (((x)>=0.0)? (x) : -(x))
  
/* work arrays size for ?levmar_der and ?levmar_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */ 
#define LM_MLE_WORKSZ(npar, nmeas) (4*(nmeas) + 4*(npar) + (nmeas)*(npar) + 2*(npar)*(npar))
  
#define LM_OPTS_SZ    	 5 
#define LM_INFO_SZ    	 10
#define LM_ERROR         -1
#define LM_INIT_MU    	 1E-03
#define LM_STOP_THRESH	 1E-20
#define LM_DIFF_DELTA    1E-20
#define LM_VERSION       "2.5 (December 2009)"
  
#define LM_CHISQ_MLE		0
#define LM_CHISQ_NEYMAN		1
#define LM_CHISQ_EQUAL_WT	2
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
		      int fitType);	/* Which type of fit to perform
					 *  LM_CHISQ_MLE is MLE for Poisson distribution
					 *  LM_CHISQ_NEYMAN is least squares where sigma^2 is set to max(x,1)
					 *  LM_CHISQ_EQUAL_WT is least squares where sigma^2 is set to 1 */
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
				      int m, int n, void *adata);
   double compute_chisq_measure (double *e, double *x, double *hx, int n,
				  int fitType, void *adata);
  char *lm_error_str (int error);
  char *fitType_str (int fitType);
  
#ifdef __cplusplus
} 
#endif				/*  */
 
#endif	/* _LEVMAR_H_ */
