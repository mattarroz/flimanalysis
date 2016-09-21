#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "fit.h"
#include "em.h"



void
iterate_amplitude_histogram (double *histogram, double ncounts,
			     Model * conv_mdl)
{
  /* comp stands for component */
  uint64_t comp, comp2, tchannel, paridx[conv_mdl->nexp];
  double sum_amp, all_comp, old_amp[conv_mdl->nexp], ncounts_inv;

  if (conv_mdl->nexp <= 1)
    {
      printf ("Internal error. The function iterate_amplitude"
	      " souldn't be called with less than or equal to"
	      " one amplitude.\n");
      return;
    }

  for (comp = 0; comp < conv_mdl->nexp; comp++)
    {
      paridx[comp] = IDX_TAU0 + 2 * comp + 1;
      old_amp[comp] = conv_mdl->par[paridx[comp]];
    }

  sum_amp = 0.0;
  ncounts_inv = 1.0 / ncounts;
  for (comp = 0; comp < conv_mdl->nexp - 1; comp++)
    {
      conv_mdl->par[paridx[comp]] = 0.0;
      for (tchannel = 0; tchannel < conv_mdl->nsamples; tchannel++)
	{
	  all_comp = 0.0;
	  for (comp2 = 0; comp2 < conv_mdl->nexp; comp2++)
	    {
	      all_comp += old_amp[comp2] * conv_mdl->pdf[comp2][tchannel];
	    }
	  /* avoid division through zero */
	  if (all_comp != 0.0)
	    {
	      conv_mdl->par[paridx[comp]] += histogram[tchannel] *
		old_amp[comp] * conv_mdl->pdf[comp][tchannel] / all_comp;
	    }
	  else
	    {
//                              printf ("Warning: division through zero.\n");
	    }
	}
      conv_mdl->par[paridx[comp]] *= ncounts_inv;
      sum_amp += conv_mdl->par[paridx[comp]];
    }

  conv_mdl->par[paridx[comp]] = 1.0 - sum_amp;

  return;
}


void
iterate_amplitude_tac (GSList * photons, uint64_t ncounts,
		       Model * conv_mdl)
{
  /* comp stands for component */
  uint64_t comp, comp2, tchannel, paridx[conv_mdl->nexp];
  double sum_amp, all_comp, old_amp[conv_mdl->nexp], ncounts_inv;
  GSList *photon;

  if (conv_mdl->nexp <= 1)
    {
      printf ("Internal error. The function iterate_amplitude"
	      " souldn't be called with less than or equal to"
	      " one amplitude.\n");
      return;
    }

  for (comp = 0; comp < conv_mdl->nexp; comp++)
    {
      paridx[comp] = IDX_TAU0 + 2 * comp + 1;
      old_amp[comp] = conv_mdl->par[paridx[comp]];
    }

  sum_amp = 0.0;

  if (!ncounts)
    {
      for (comp = 0; comp < conv_mdl->nexp - 1; comp++)
	{
	  conv_mdl->par[paridx[comp]] = 0.0;
	}
      conv_mdl->par[paridx[comp]] = 1.0;

      return;
    }
  ncounts_inv = 1.0 / ncounts;

  for (comp = 0; comp < conv_mdl->nexp - 1; comp++)
    {
      conv_mdl->par[paridx[comp]] = 0.0;
      for (photon = photons; photon; photon = g_slist_next (photon))
	{
	  all_comp = 0.0;
	  tchannel = GPOINTER_TO_UINT (photon->data);
	  //printf ("%d ", tchannel);

	  for (comp2 = 0; comp2 < conv_mdl->nexp; comp2++)
	    {
	      all_comp += old_amp[comp2] * conv_mdl->pdf[comp2][tchannel];
	    }
	  /* avoid division through zero */
	  if (all_comp != 0.0)
	    {
	      conv_mdl->par[paridx[comp]] +=
		old_amp[comp] * conv_mdl->pdf[comp][tchannel] / all_comp;
	    }
	  else
	    {
//                              printf ("Warning: division through zero.\n");
	    }
	}
      conv_mdl->par[paridx[comp]] *= ncounts_inv;
      sum_amp += conv_mdl->par[paridx[comp]];
    }

  conv_mdl->par[paridx[comp]] = 1.0 - sum_amp;

  return;
}
