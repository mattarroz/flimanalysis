#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics.h>

#include "gnuplot_i.h"
#include "utils.h"

double
norm_histogram (double *histogram, uint64_t ntchannels, double ncounts)
{
  uint64_t tchannel;
  double integral, norm_factor;

  integral = 0.0;
  for (tchannel = 0; tchannel < ntchannels; tchannel++)
    {
      integral += histogram[tchannel];
    }

  norm_factor = ncounts / integral;
  for (tchannel = 0; tchannel < ntchannels; tchannel++)
    {
      histogram[tchannel] *= norm_factor;
    }

  return norm_factor;
}

double
norm_histogram_peak (double *histogram, uint64_t ntchannels, uint64_t ncounts)
{
  uint64_t tchannel;
  double max, factor;

  max = 0.0;
  for (tchannel = 0; tchannel < ntchannels; tchannel++)
    {
      if (histogram[tchannel] > max)
        max = histogram[tchannel];
    }

  factor = ncounts / max;
  for (tchannel = 0; tchannel < ntchannels; tchannel++)
    {
      histogram[tchannel] *= factor;
    }

  return factor;
}

int
shift_histogram (double *a, int N, double s)
{

  double mant = 0.0;
  int i = 0, si = 0;
  if (!s || (s != s))
    return 1;
  if ((s > (N - 1)) || (s < -(N - 1)))
    return 0;

  mant = s - floor (s);
  si = (int) (s - mant);

  if (si < 0)
    {
      si = -si;
      for (i = 0; i < (N - si); i++)
        a[i] = (a[i + si - 1] - a[i + si]) * mant + a[i + si];
      //    assert(((i+si)<=N));
      a[N - si] = a[N - 1] * mant;
      for (i = N - si + 1; i < N; i++)
        a[i] = 0.0;
      //    assert((i<=N));
    }
  else
    {
      if (si > 0)
        {
          for (i = N - 1; i > si; i--)
            a[i] = (a[i - si - 1] - a[i - si]) * mant + a[i - si];
          //    assert(i-si-1>=-1);
          a[si] = a[0] * mant;
          for (i = 0; i < (si); i++)
            a[i] = 0.0;
        }
      else
        {
          for (i = N - 1; i >= 1; i--)
            a[i] = (a[i - 1] - a[i]) * mant + a[i];
          //       assert(i<N+1);
          a[0] = a[0] * (1 - mant);
        }

    }
//if(SHIDEBUG){for(i=0;i<N;i++)
  //       printf("%f ",a[i]);}
  return 1;
}

double
calc_ncounts (double *histogram, uint64_t size)
{
  double ncounts;
  uint64_t i;

  ncounts = 0;
  for (i = 0; i < size; i++)
    ncounts += histogram[i];

  return ncounts;
}

void
plot_double (double *sample, uint64_t ntchannels)
{
  gnuplot_ctrl *ploth;

  ploth = gnuplot_init ();
  gnuplot_plot_x (ploth, sample, ntchannels, "data");
  getchar ();
  gnuplot_close (ploth);

  return;
}

void
plot_double_tex (double *sample, char *xtitle, char *ytitle, uint64_t ntchannels)
{
  gnuplot_ctrl *ploth;
  char fname[255] = "out", buf[255];
  double border_space;

  ploth = gnuplot_init ();
  sprintf (buf, "\"%s.tex\"", fname);
  setup_gnuplot_tex (ploth, buf, NULL);
  gnuplot_set_xlabel (ploth, xtitle);
  gnuplot_set_ylabel (ploth, ytitle);

  border_space = ntchannels * 0.1;
  sprintf (buf, "set xrange [%f:%f]", 0 - border_space,
           ntchannels + border_space);
  gnuplot_cmd (ploth, buf);
  gnuplot_plot_x (ploth, sample, ntchannels, "data");

  gnuplot_close (ploth);

  sprintf (buf, "latex %s", fname);
  system (buf);
  sprintf (buf, "dvipdf %s.dvi", fname);
  system (buf);
  sprintf (buf, "rm -f %s.tex %s.dvi %s.log %s-inc.eps %s.aux q.log", fname,
           fname, fname, fname, fname);
  system (buf);

  return;
}

void
plot_double_xy_tex (double *x, char *xtitle, double *y, char *ytitle,
                    uint64_t ntchannels)
{
  gnuplot_ctrl *ploth;
  char fname[255] = "out", buf[255];
  double border_space;

  ploth = gnuplot_init ();
  sprintf (buf, "\"%s.tex\"", fname);
  setup_gnuplot_tex (ploth, buf, NULL);
  gnuplot_set_xlabel (ploth, xtitle);
  gnuplot_set_ylabel (ploth, ytitle);

  border_space =
      fabs (gsl_stats_max (x, 1, ntchannels) -
            gsl_stats_min (x, 1, ntchannels)) * 0.1;
  sprintf (buf, "set xrange [%f:%f]",
           gsl_stats_min (x, 1, ntchannels) - border_space, gsl_stats_max (x,
                                                                           1,
                                                                           ntchannels)
                                                            + border_space);
  gnuplot_cmd (ploth, buf);

  border_space =
      fabs (gsl_stats_max (y, 1, ntchannels) -
            gsl_stats_min (y, 1, ntchannels)) * 0.1;
  sprintf (buf, "set yrange [%f:%f]",
           gsl_stats_min (y, 1, ntchannels) - border_space, gsl_stats_max (y,
                                                                           1,
                                                                           ntchannels)
                                                            + border_space);
  gnuplot_cmd (ploth, buf);

  gnuplot_plot_xy (ploth, x, y, ntchannels, "data");

  gnuplot_close (ploth);

  sprintf (buf, "latex %s", fname);
  system (buf);
  sprintf (buf, "dvipdf %s.dvi", fname);
  system (buf);
  sprintf (buf, "rm -f %s.tex %s.dvi %s.log %s-inc.eps %s.aux q.log", fname,
           fname, fname, fname, fname);
  system (buf);

  return;
}

void
plot_double_xy_yerror_tex (double *x, char *xtitle, double *y, char *ytitle,
                           double *yerror, uint64_t ntchannels)
{
  gnuplot_ctrl *ploth;
  char fname[255] = "out", buf[255];
  uint64_t i;
  FILE *pFile;
  FILE *gnuFile = fopen ("./plot.gnu", "w");

  ploth = gnuplot_init ();
  sprintf (buf, "\"%s.tex\"", fname);
  setup_gnuplot_tex (ploth, buf, gnuFile);
  gnuplot_set_xlabel (ploth, xtitle);
  fprintf (gnuFile, "set xtitle '%s'\n", xtitle);
  gnuplot_set_ylabel (ploth, ytitle);
  fprintf (gnuFile, "set ytitle '%s'\n", ytitle);

  gnuplot_cmd (ploth, "set logscale x");
  fprintf (gnuFile, "set logscale x");

  pFile = fopen ("/tmp/temp.dat", "w");
  for (i = 0; i < ntchannels; i++)
    {
      fprintf (pFile, "%f %f %f\n", x[i], y[i], yerror[i]);
    }
  fclose (pFile);
  gnuplot_cmd (ploth, "plot \"/tmp/temp.dat\" with yerrorbars");
  fprintf (gnuFile, "plot \"/tmp/temp.dat\" with yerrorbars");

  gnuplot_close (ploth);
  fclose (gnuFile);

  sprintf (buf, "latex %s", fname);
  system (buf);
  sprintf (buf, "dvipdf %s.dvi", fname);
  system (buf);
  sprintf (buf, "rm -f %s.tex %s.dvi %s.log %s-inc.eps %s.aux q.log", fname,
           fname, fname, fname, fname);
  system (buf);

  return;
}

void
setup_gnuplot_tex (gnuplot_ctrl *ploth, char *fname, FILE *gnuFile)
{
  char *set_output_cmd, *set_output = "set output ";

  gnuplot_cmd (ploth,
               "set term epslatex standalone color header \" \\\\usepackage[T1]{fontenc} \\\\usepackage{fourier} \\\\usepackage{amsmath,amsfonts,amssymb} \"");
  set_output_cmd =
      malloc (sizeof (char) * (strlen (set_output) + strlen (fname) + 1));
  sprintf (set_output_cmd, "%s%s", set_output, fname);
  gnuplot_cmd (ploth, set_output_cmd);
  gnuplot_cmd (ploth, "unset key");
  free (set_output_cmd);

  if (gnuFile != NULL)
    {
      fprintf (gnuFile,
               "set term epslatex standalone color header \" \\\\usepackage[T1]{fontenc} \\\\usepackage{fourier} \\\\usepackage{amsmath,amsfonts,amssymb} \"");
      fprintf (gnuFile, "%s\n", set_output_cmd);
      fprintf (gnuFile, "unset key\n");
    }
}
