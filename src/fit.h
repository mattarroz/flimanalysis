#include <stdint.h>
#include <glib.h>

typedef struct _Histogram
{
  double ncounts;
  double ncounts_peak;
  uint16_t nsamples;
  double tcal;
  double *irf;
  double *sample;
  double *darknoise;
} Histogram;

typedef struct _Pixel
{
  uint16_t x;
  uint16_t y;
  GSList *photons;
  uint32_t ncounts;
} Pixel;


// Not used at the moment
typedef struct _Photon
{
  uint16_t TAC;
  uint64_t T;

} Photon;


/* Structure containing the data which is handed over to the convolution
	algorithm */
typedef struct _Model
{
  uint8_t nexp;
  uint16_t tstart;
  uint16_t tstop;

  double ncounts;
  double norm_factor;
#define IDX_BG		0
#define IDX_IRFSHIFT	1
#define IDX_TAU0	2
  double *par;
//    double **link_par;
  uint16_t npar;

#define PAR_VARIABLE	-2
#define PAR_FIXED	-1
  uint16_t *pvar;
  uint16_t npvar;

  uint16_t nsamples;		// number of time channels

  // TODO: alle referenzen auf irf sollten in irf_roi reinwandern
  double *irf;
  double *shifted_irf;
  double shifted_irf_offset;

  uint16_t irf_tstart;
  uint16_t irf_tstop;


  double *darknoise;
#define YES		1
#define NO		0
  uint16_t constant_darknoise;	/* if set to YES, we have constant background
				   over all time channels and the dark
				   measurement is unused */

  uint8_t fitType;

  double **pdf;			/* single 1-normalized pdfs, for all 
				   lifetimes and pixels */
  double *complete_pdf;		/* complete 1-normalized pdf, one for 
				   every pixel */
  double *fit_function;		/* fit function (completed pdf normalized
				   to the number of counts), one for every 
				   pixel */
} Model;


typedef struct _ROI
{
  uint64_t npixels;
  Pixel *pixel;
  double *histogram;
  Model *mdl;
  uint64_t ncounts;

  // TODO: entfernen
  uint16_t width;
  uint16_t height;

  uint16_t tstart;
  uint16_t tstop;
} ROI;


typedef struct _Reference
{
  uint16_t iroi;			/* Index of ROI */
  uint16_t ipar;			/* Index of fixed parameter */
  uint16_t ipvar;			/* Index of variable parameter */
  uint16_t ilink;			/* set to zero if parameter not linked */
} Reference;


typedef struct _Global_Model
{
  uint8_t nexp;
  uint16_t nlinks;
  double *link;

  /* pvar_multi_map and pvar_multi share the same index */
  GList *pvar_multi_map;
  double *pvar_multi;
  uint16_t npvar;
  uint8_t fitType;
  uint16_t threshold;
} Global_Model;

typedef struct _FLIM_measurement
{
  double *joint_histogram;	// histograms of all ROIs
  ROI *roi;
  uint16_t nrois;
  GSList ***pixellist;
  GSList ***pixelindex;
  uint16_t nsamples;
  uint16_t width;
  uint16_t height;
  Global_Model *global_mdl;
  uint64_t ncounts;

  ROI *irf_roi;
  uint16_t nirf_rois;
  uint64_t ncounts_irf;
  GSList ***irf_pixellist;
  GSList ***irf_pixelindex;
  double *joint_histogram_irf;

  double tcal;
  double tabs_clock;
  uint16_t pixeldwelltime;
  uint64_t tabs_start;
  uint64_t tabs_stop;

  char *meas_fname;
  char *irf_fname;
  uint8_t fileformat;
#define QADATA 1
#define PT3DATA 2
} FLIM_measurement;



uint8_t fit_levmar_FLIM_measurement (FLIM_measurement * m);
uint8_t fit_levmar_histogram (double *histogram, Model * mdl, uint8_t fitType);
uint8_t fit_exp_gs (Histogram * measure, Model * mdl, uint8_t fitType);
void fit_FLIM_measurement (FLIM_measurement * m);
double *cramer_rao_boundary_gs (Model * mdl);

// return values for em algorithm
#define	SMALL_DELTA	0
#define MAX_IT		1
uint16_t fit_amplitude_histogram (double *histogram, double ncounts,
			      Model * conv_mdl);
uint16_t fit_amplitude_tac (GSList * photons, uint64_t ncounts,
			Model * conv_mdl);
char *emerr2str (uint8_t res);

ROI *create_rois (uint16_t nrois, uint16_t nsamples, uint8_t nexp);

void free_measurement (Histogram * measure);	// TODO: in free_histogram umbennenen
void free_convolution_par (Model * conv_mdl);
void free_FLIM_measurement (FLIM_measurement * m);
void free_ROI (ROI * roi);
void free_Pixel (Pixel * pixel);
