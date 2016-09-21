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
  
write_file (double *p, int n) 
{
  
  
  
					   a file or create a file if it does not exist. */
  
    
  
  



compute_chisq_measure (double *e, double *x, double *hx, int n, int fitType,
		       void *adata) 
{
  
  
  

    
    {
    
      
	
	{
	  
	    
	    {
	      
		
	      
	      else
		
	      
	    
	  
	  else if (hx[i] != hx[i])
	    {
	      
	    
	  
	  else
	    
	
      
    
      
	
	{
	  
	    {
	      
		
		{
		  
		  
		
	      
	      else
		
		{
		  
		  
		
	    
	  else
	    {
	      
	    
	
      
    
      
	
	{
	  
	  
	
      
    
  



dlevmar_mle_fdif_forw_jac_approx (
				  (double *p, double *hx, int m, int n,
				   void *adata), 
				  /* function to differentiate */ 
				  double *p, /* I: current parameter estimate, mx1 */ 
				  double *hx, /* I: func evaluated at p, i.e. hx=func(p), nx1 */ 
				  double *hxx, /* W/O: work array for evaluating func(p+delta), nx1 */ 
				  double delta, /* increment for computing the Jacobian */ 
				  double *jac, /* O: array for storing approximated Jacobian, nxm */ 
				  int m, 
{
  
  
  
  
//double dif[4096];
    
    {
      
	/* determine d=max(1E-04*|p[j]|, delta), see HZ */ 
	d = 1E-4 * p[j];	// force evaluation
      d = FABS (d);
      
	
      
      
      
      
      
      
	{
	  
	  
//      dif[i] = hxx[i]-hx[i];
	}
      
//    plot_double (jac, n*m);
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

  dlevmar_mle_diff (
		    (double *p, double *hx, int m, int n, void *adata),
		    /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
		    
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
  
    tmp;			/* mainly used in matrix & vector multiplications */
  
  
  
  
  
  
  const int nm = n * m;
  
  
  
  
  
  
    {
      
		"(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n",
		n, m);
      
    
  
    //if(!jacf){
    //  fprintf(stderr, "No function specified for computing the Jacobian");
    //  return LM_ERROR;
    //}
    
    {
      
      
      
      
      
      
    
  
  else
    {				// use default values
      tau = LM_INIT_MU;
      
      
      
      
      
    
  
    {
      
      work = (double *) malloc (worksz * sizeof (double));	/* allocate a big chunk in one step */
      
	{
	  
	  
	
      
    
  
    /* set up work arrays */ 
    e = work;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    /* compute e=x - f(p) and its chisq measure */ 
    (*func) (p, hx, m, n, adata);
  nfev = 1;
  
    /* ### e=x-hx, p_eL2=||e|| */ 
    p_chisq = compute_chisq_measure (e, x, hx, n, fitType, adata);
  
  
    stop = 7;
  
    {
      
      
      
	/* Note that p and e have been updated at a previous iteration */ 
	
	{			/* error is small */
	  
	  
	
      
	/* Compute the Jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
	 * Since J^T J is symmetric, its computation can be sped up by computing
	 * only its upper triangular part and copying it to the lower part
	 */ 
	
	//(*jacf)(p, jac, m, n, adata); ++njev;
	//if((updp && nu>16) || updjac==K){ /* compute difference approximation to J */
	
					     n, adata);
      
	//write_file(jac,n);
//      plot_double (jac, n);
	nfev += m;
      
	//nu=2;
	//}
	
	/* J^T J, J^T e */ 
	for (i = m * m; i-- > 0;)
	
      
	
      
	{
	  
	  
	    
	    {
	    
	      
		
	      
	      else
		
	      
	    
	      
		
	      
	      else
		
	      
	    
	      
	      
	    
	  
	    {
	      
	      
	      for (j = i + 1; j-- > 0;)	/* j<=i computes lower triangular part only */
		
	      
		/* J^T e */ 
		jacTe[i] += alpha * e[l];
	    
	
      
	
	  
      
	/* Compute ||J^T e||_inf and ||p||^2 */ 
	for (i = 0, p_L2 = jacTe_inf = 0.0; i < m; ++i)
	{
	  
	    jacTe_inf = tmp;
	  
	  
	
      
	/* check for convergence */ 
	if ((jacTe_inf <= eps1))
	{
	  
	  
	  
	
      
	/* compute initial damping factor */ 
	if (k == 0)
	{
	  
	    
	      tmp = diag_jacTjac[i];	/* find max diagonal element */
	  
	
      
	/* determine increment using adaptive damping */ 
	while (1)
	{
	  
	    /* augment normal equations */ 
	    for (i = 0; i < m; ++i)
	    
	  
	    /* solve augmented equations */ 
	    /* use the LU included with GSL */ 
	    for (i = 0; i < m * m; i++)
	    
	  
	  
	    !gsl_linalg_LU_solve (&LUdecomp_view.matrix, perm,
				  &jacTe_view.vector, &Dp_view.vector);
	  
	    {
	      
		/* compute p's new estimate and ||Dp||^2 */ 
		for (i = 0, Dp_L2 = 0.0; i < m; ++i)
		{
		  
		  
		
	      
		//Dp_L2=sqrt(Dp_L2);
		
		{		/* relative change in p is small, stop */
		  
		  
		
	      
//       if(Dp_L2>=(p_L2+eps2)/EPSILON*EPSILON){ /* almost singular */
//         stop=4;
//         break;
//      }
		
	      ++nfev;		/* evaluate function at p + Dp */
	      
		/* compute ||e(pDp)||_2 */ 
		/* ### hx=x-hx, pDp_chisq=||hx|| */ 
		pDp_chisq =
		compute_chisq_measure (e_test, x, hx, n, fitType, adata);
	      
		{		/* chisq is not finite, most probably due to a user error.
				 * This check makes sure that the inner loop does not run indefinitely.
				 * Thanks to Steve Danauskas for reporting such cases
				 */
		  
		  
		
	      
		
	      
	      
		{		/* reduction in error, increment is accepted */
		  
		  
		  
		  
		  
		    
		  
		    
		  
		  
		    
		  
		
	    
	  
	    /* if this point is reached, either the linear system could not be solved or
	     * the error did not reduce; in any case, the increment must be rejected
	     */ 
	    
	  
	  if (nu2 <= nu)
	    {			/* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
	      
	      
	    
	  
	  
	    
	
    
  
    stop = 3;
  
    
  
    {
      
      
      
      
      
	
	  tmp = jacTjac[i * m + i];
      
      
      
      
      
      
    
  
    free (work);
  
  



lm_error_str (int error)
{
  
  
  
    {
    
      
      
    
      
      
    
      
      
    
      
	"singular matrix. Restart from current p with increased mu";
      
    
      
	"no further error reduction is possible. Restart with increased mu";
      
    
      
      
    
      
	"stopped by invalid (i.e. NaN or Inf) \"func\" values. This is a user error";
      
    
      
	"stopped by small change in chisq on successful iteration (dF<eps4)";
      
    
      
      
    
  



fitType_str (int fitType)
{
  
  
    {
    
      
      
    
      
      
    
      
      
    
      
    
  


