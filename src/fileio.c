#include <stdio.h>
#include <math.h>
#include <tiffio.h>

#include "fit.h"
#include "fileio.h"
#include "utils.h"
#include "ini.h"
#include "readpt3.h"
#include "read_qadata.h"
#include "levmar_mle_diff.h"
#include "convolution.h"
//#include "roireader_free/ROILib.h"

//#define DEBUG 1


// TODO: free !!!! --> sollte ich mit einer "Aufräumfunktion" lösen, die auch bei allen fehlern dann benutzt wird


int checkfile (char file[], uint16_t * rowstart, uint16_t * rowend, double *timestep);
void readfile (char file[], uint64_t rowstart, uint64_t length,
		double *fluorescence, double *irf);
int calculate_roi_pixels (uint16_t xbottom, uint16_t ybottom, uint16_t xtop, uint16_t ytop,
		ROI *roi, FLIM_measurement * m);
void calculate_roi_pixels_from_roi_file (char *filename, uint64_t raw_roi_index,
		uint16_t width, uint16_t height,
		uint16_t xbottom, uint16_t ybottom,
		uint16_t xtop, uint16_t ytop, uint16_t iroi,
		FLIM_measurement * m);
void create_pixellist (FLIM_measurement * m);
void create_pixellist_from_TIFF (char *filename, FLIM_measurement * m);
void set_link_fixed (FLIM_measurement * m, uint16_t ipar, char *keyvalue_raw,
		uint16_t iroi, uint8_t islifetime);
void set_variable (FLIM_measurement * m, uint16_t ipar, char *keyvalue_raw,
		uint16_t iroi);

/* TODO: 
	- fehlerüberprüfung? ich kann ins ini file werte eingeben, die das programm schrotten FIXME!!!!!

 */

FLIM_measurement *
read_flim_measurement (char *ini_file)
{
	uint64_t res;
	struct CfgStruct cfg;
	char roifile_tiff[255]="";
	FLIM_measurement *m;
	char sectionname[255], buf[255];
	uint16_t xtop, ytop, ybottom, xbottom, i, j;
	FILE *CfgFile;
	char keyname[255], keyvalue[255];

	m = malloc (sizeof (FLIM_measurement));

	if (NULL == (CfgFile = fopen (ini_file, "r")))
	{
		printf ("Error. Could not open file %s.\n", ini_file);
		return NULL;
	}
	// Section [Globals]
	cfg.Name = "dataformat";
	cfg.DataPtr = &keyvalue;
	cfg.VarType = Cfg_String;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	if (strcmp (keyvalue, "QA") == 0) {
		m->meas_fname = calloc (256,sizeof(char));
	    cfg.Name = "qadata";
	    cfg.DataPtr = m->meas_fname;
	    cfg.VarType = Cfg_String;
	    ReadCfg(ini_file, "Globals", &cfg, CfgFile);

	    m->fileformat = QADATA;
	} else if (strcmp (keyvalue, "pt3") == 0) {
		m->meas_fname = calloc (256,sizeof(char));
		cfg.Name = "pt3data";
		cfg.DataPtr = m->meas_fname;
		cfg.VarType = Cfg_String;
		ReadCfg (ini_file, "Globals", &cfg, CfgFile);

		m->fileformat = PT3DATA;
	} else {
		printf ("Error: key \"%s\" was not found or equal \"pt3\" or \"QA\" (upper case) in file %s. \n.", cfg.Name,
				ini_file);
		return NULL;
	}
	m->irf_fname = calloc (256,sizeof(char));
	cfg.Name = "irf";
	cfg.DataPtr = m->irf_fname;
	cfg.VarType = Cfg_String;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "width";
	cfg.DataPtr = &(m->width);
	cfg.VarType = Cfg_Uint16;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "height";
	cfg.DataPtr = &(m->height);
	cfg.VarType = Cfg_Uint16;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "samples";
	cfg.DataPtr = &(m->nsamples);
	cfg.VarType = Cfg_Uint16;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "rois";
	cfg.VarType = Cfg_String;
	cfg.DataPtr = buf;
	if (ReadCfg (ini_file, "Globals", &cfg, CfgFile) < 1)
	{
		printf ("Error while reading \"rois\" key in %s\n.", ini_file);
		return NULL;
	}
	if (strcmp (buf, "file") == 0)
	{
		cfg.Name = "roifile_tiff";
		cfg.DataPtr = roifile_tiff;
		if (ReadCfg (ini_file, "Globals", &cfg, CfgFile) < 1)
		{
			printf ("Error while reading \"%s\" key in %s\n.", cfg.Name,
					ini_file);
			return NULL;
		}
	} else {
		cfg.DataPtr = &(m->nrois);
		cfg.VarType = Cfg_Uint16;
		res = ReadCfg (ini_file, "Globals", &cfg, CfgFile);
		if (m->nrois == 0) {
				printf ("At least one ROI needs to be defined.\n");
				return NULL;
			}
	}


	cfg.Name = "irfrois";
	cfg.DataPtr = &(m->nirf_rois);
	cfg.VarType = Cfg_Uint16;
	if (ReadCfg (ini_file, "Globals", &cfg, CfgFile) < 1)
	{
		printf ("Error while reading \"irfrois\" key in %s\n.", ini_file);
		return NULL;
	}
	if (m->nirf_rois == 0) {
			printf ("At least one ROI needs to be defined.\n");
			return NULL;
		}
	m->joint_histogram_irf =
			calloc (m->nsamples * m->nirf_rois, sizeof (double));


	if (m->fileformat == QADATA) {
		m->irf_roi = calloc(m->nirf_rois, sizeof(ROI));
		for (i = 0; i < m->nirf_rois; i++) {
			m->irf_roi[i].histogram =
					&(m->joint_histogram_irf[i * m->nsamples]);

			sprintf(sectionname, "IRFROI%d", i);

			cfg.Name = "tstart";
			cfg.DataPtr = &m->irf_roi[i].tstart;
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "tstop";
			cfg.DataPtr = &m->irf_roi[i].tstop;
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "rectangle.xtop";
			cfg.DataPtr = &(xtop);
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "rectangle.ytop";
			cfg.DataPtr = &(ytop);
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "rectangle.xbottom";
			cfg.DataPtr = &(xbottom);
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "rectangle.ybottom";
			cfg.DataPtr = &(ybottom);
			cfg.VarType = Cfg_Uint16;
			ReadCfg(ini_file, sectionname, &cfg, CfgFile);

			if (!calculate_roi_pixels(xbottom, ybottom, xtop, ytop, &(m->irf_roi[i]), m)) {
				// TODO: free everything
				return NULL;
			}
#ifdef DEBUG
			printf("xtop: %d xbottom: %d ytop: %d ybottom: %d\n", xtop,
					xbottom, ytop, ybottom);
#endif
		}
	}

	cfg.Name = "tcal";
	cfg.DataPtr = &(m->tcal);
	cfg.VarType = Cfg_Double;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "tabsclock";
	cfg.DataPtr = &(m->tabs_clock);
	cfg.VarType = Cfg_Double;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);
	m->tabs_clock *= 1.0e3;

	cfg.Name = "tabsstart";
	cfg.DataPtr = &(m->tabs_start);
	cfg.VarType = Cfg_Uint64;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "tabsstop";
	cfg.DataPtr = &(m->tabs_stop);
	cfg.VarType = Cfg_Uint64;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	m->global_mdl = malloc (sizeof (Model));
	cfg.Name = "nexp";
	cfg.DataPtr = &(m->global_mdl->nexp);
	cfg.VarType = Cfg_Uint8;
	ReadCfg (ini_file, "Model", &cfg, CfgFile);
#ifdef DEBUG
	printf ("nexp %d\n", m->global_mdl->nexp);
#endif

	cfg.Name = "threshold";
	cfg.DataPtr = &(m->global_mdl->threshold);
	cfg.VarType = Cfg_Uint16;
	ReadCfg (ini_file, "Globals", &cfg, CfgFile);

	cfg.Name = "pixeldwelltime";
	cfg.DataPtr = &(m->pixeldwelltime);
	cfg.VarType = Cfg_Uint16;
	if (ReadCfg (ini_file,"Globals", &cfg, CfgFile) < 1) {
		printf ("Could not read pixel dwell time from %s\n", ini_file);
		return NULL;
	}


	cfg.Name = "nlinks";
	cfg.DataPtr = &(m->global_mdl->nlinks);
	cfg.VarType = Cfg_Uint16;
	ReadCfg (ini_file, "Model", &cfg, CfgFile);
#ifdef DEBUG
	printf ("nlinks %d\n", m->global_mdl->nlinks);
#endif

	m->global_mdl->link = calloc (m->global_mdl->nlinks, sizeof (double));
	for (i = 0; i < m->global_mdl->nlinks; i++)
	{
		sprintf (keyname, "link[%d]", i);

		cfg.Name = keyname;
		cfg.DataPtr = &(m->global_mdl->link[i]);
		cfg.VarType = Cfg_Double;
		ReadCfg (ini_file, "Model", &cfg, CfgFile);

#ifdef DEBUG
		printf ("link[%d]: %f\n", i, m->global_mdl->link[i]);
#endif
	}


	cfg.Name = "fittype";
	cfg.DataPtr = &buf;
	cfg.VarType = Cfg_String;
	ReadCfg (ini_file, "Model", &cfg, CfgFile);

	if (strcmp (buf, "mle") == 0)
	{
		m->global_mdl->fitType = LM_CHISQ_MLE;
	}
	else if (strcmp (buf, "chisq_equal_wt") == 0)
	{
		m->global_mdl->fitType = LM_CHISQ_EQUAL_WT;
	}
	else if (strcmp (buf, "chisq_neyman") == 0)
	{
		m->global_mdl->fitType = LM_CHISQ_NEYMAN;
	}
	else
	{
#ifdef DEBUG
		printf ("Could not read key fitType. Setting fitType to MLE\n");
#endif
		m->global_mdl->fitType = LM_CHISQ_MLE;
	}

	if (strcmp(roifile_tiff,"") != 0)
	{

		if (!calculate_roi_pixels_from_tiff_file(roifile_tiff,m)) {
			return NULL;
		}

		m->global_mdl->pvar_multi =
			calloc ((IDX_TAU0 + m->global_mdl->nexp + 1) * m->nrois, sizeof (double));

		m->global_mdl->pvar_multi_map = NULL;
		m->global_mdl->npvar = 0;


		// soll ich überhaupt unterschiedliche fitparameter für unterschiedliche rois erlauben?
		// bei den TIFFs die macht das eigentlich keinen Sinn, wenn man eh nicht weiß
		// welche ROI zu was genau gehört...
		sprintf (sectionname, "Model");

		for (i = 0; i < m->nrois; i++) {

			// We don't allow different fitting algorithms for different ROIs
			m->roi[i].mdl->fitType = m->global_mdl->fitType;
			m->roi[i].mdl->irf = m->joint_histogram_irf;
			m->roi[i].histogram = &(m->joint_histogram[i * m->nsamples]);

			cfg.Name = "tstart";
			cfg.DataPtr = &(m->roi[i].mdl->tstart);
			cfg.VarType = Cfg_Uint16;
			ReadCfg (ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "tstop";
			cfg.DataPtr = &(m->roi[i].mdl->tstop);
			cfg.VarType = Cfg_Uint16;
			ReadCfg (ini_file, sectionname, &cfg, CfgFile);

			// IRFSHIFT + BACKGROUND

			strcpy (keyname, "irfshift");
			cfg.Name = keyname;
			cfg.DataPtr = &(m->roi[i].mdl->par[IDX_IRFSHIFT]);
			cfg.VarType = Cfg_Double;
			res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
			sprintf (keyname, "link[%d]",m->global_mdl->nexp);
	#ifdef DEBUG
			printf ("%s: %f\n", keyname, m->roi[i].mdl->par[IDX_IRFSHIFT]);
	#endif
			set_link_fixed (m, IDX_IRFSHIFT, keyname, i, 0);

			strcpy (keyname, "background");
			cfg.Name = keyname;
			cfg.DataPtr = &(m->roi[i].mdl->par[IDX_BG]);
			cfg.VarType = Cfg_Double;
			res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
	#ifdef DEBUG
			printf ("%s: %f\n", keyname, m->roi[i].mdl->par[IDX_BG]);
	#endif
			set_link_fixed (m, IDX_BG, keyname, i, 0);

			for (j = 0; j < m->global_mdl->nexp; j++)
			{


				sprintf (keyname, "lifetime[%d]", j);
				cfg.Name = keyname;
				cfg.DataPtr = &(m->roi[i].mdl->par[IDX_TAU0 + 2 * j]);
				cfg.VarType = Cfg_Double;
				res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
				m->roi[i].mdl->par[IDX_TAU0 + 2 * j] *= 1.0 / m->tcal;
				sprintf (keyname, "link[%d]",j);
#ifdef DEBUG
				printf ("%s: %f\n", keyname,
						m->roi[i].mdl->par[IDX_TAU0 + 2 * j]);
#endif
				set_link_fixed (m, IDX_TAU0 + 2 * j, keyname, i, 1);

				sprintf (keyname, "amp[%d]", j);
				cfg.Name = keyname;
				cfg.DataPtr = &(m->roi[i].mdl->par[IDX_TAU0 + 2 * j + 1]);
				cfg.VarType = Cfg_Double;
				res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
#ifdef DEBUG
				printf ("%s: %f\n", keyname,
						m->roi[i].mdl->par[IDX_TAU0 + 2 * j + 1]);
#endif
				if (res < 1)
				{
					cfg.DataPtr = keyvalue;
					cfg.VarType = Cfg_String;
					ReadCfg (ini_file, sectionname, &cfg, CfgFile);
					set_link_fixed (m, IDX_TAU0 + 2 * j + 1, keyvalue, i, 0);
				}
				else if (((j + 1) == m->global_mdl->nexp)
						&& (m->global_mdl->fitType == LM_CHISQ_MLE))
				{
#ifdef DEBUG
					printf
					("last amplitude par[%d] given through normalization\n",
							IDX_TAU0 + 2 * j + 1);
#endif
				}
				else
				{
					set_variable (m, IDX_TAU0 + 2 * j + 1, keyvalue, i);
				}
			}

		}
	}
	else
	{
		m->global_mdl->pvar_multi =
			calloc ((IDX_TAU0 + m->global_mdl->nexp + 1) * m->nrois, sizeof (double));

		m->global_mdl->pvar_multi_map = NULL;
		m->global_mdl->npvar = 0;
		m->roi = create_rois (m->nrois, m->nsamples, m->global_mdl->nexp);

		m->joint_histogram = calloc (m->nrois * m->nsamples, sizeof (double));
#ifdef DEBUG
		printf ("nrois %d\n", m->nrois);
#endif

		for (i = 0; i < m->nrois; i++)
		{
			char keyname[255], keyvalue[255];
			int res;

			m->roi[i].mdl->fitType = m->global_mdl->fitType;
			m->roi[i].histogram = &(m->joint_histogram[i * m->nsamples]);
			m->roi[i].mdl->irf = m->joint_histogram_irf;
			m->roi[i].histogram = &(m->joint_histogram[i * m->nsamples]);

			sprintf (sectionname, "Model");
			cfg.Name = "tstart";
			cfg.DataPtr = &(m->roi[i].mdl->tstart);
			cfg.VarType = Cfg_Uint16;
			ReadCfg (ini_file, sectionname, &cfg, CfgFile);

			cfg.Name = "tstop";
			cfg.DataPtr = &(m->roi[i].mdl->tstop);
			cfg.VarType = Cfg_Uint16;
			ReadCfg (ini_file, sectionname, &cfg, CfgFile);

			sprintf (sectionname, "ROI%d", i);

			cfg.Name = "rectangle.xtop";
			cfg.DataPtr = &(xtop);
			cfg.VarType = Cfg_Uint16;
			if (ReadCfg (ini_file, sectionname, &cfg, CfgFile) < 1)
			{
				xtop = 0;
			}

			cfg.Name = "rectangle.ytop";
			cfg.DataPtr = &(ytop);
			cfg.VarType = Cfg_Uint16;
			if (ReadCfg (ini_file, sectionname, &cfg, CfgFile) < 1)
			{
				ytop = 0;
			}

			cfg.Name = "rectangle.xbottom";
			cfg.DataPtr = &(xbottom);
			cfg.VarType = Cfg_Uint16;
			if (ReadCfg (ini_file, sectionname, &cfg, CfgFile) < 1)
			{
				xbottom = m->width - 1;
			}

			cfg.Name = "rectangle.ybottom";
			cfg.DataPtr = &(ybottom);
			cfg.VarType = Cfg_Uint16;
			if (ReadCfg (ini_file, sectionname, &cfg, CfgFile) < 1)
			{
				ybottom = m->height - 1;
			}



#ifdef DEBUG
			printf ("xtop: %d xbottom: %d ytop: %d ybottom: %d\n", xtop,
					xbottom, ytop, ybottom);
#endif

			m->roi[i].width = xbottom - xtop + 1;
			m->roi[i].height = ybottom - ytop + 1;

			cfg.Name = "roitype";
			cfg.DataPtr = &(keyvalue);
			cfg.VarType = Cfg_String;
			ReadCfg (ini_file, sectionname, &cfg, CfgFile);
			if (strcmp (keyvalue, "rectangle") == 0)
			{
				if (!calculate_roi_pixels (xbottom, ybottom, xtop, ytop, &(m->roi[i]), m))
				{
					printf ("Failed while calculating roi pixels for roi %d\n", i);
					return NULL;
				}
			}
			else
			{
				printf ("Error while reading roitype of ROI %d. Exiting...\n",
						i);
				return NULL;
			}

			for (j = 0; j < m->global_mdl->nexp; j++)
			{
				sprintf (keyname, "lifetime[%d]", j);
				cfg.Name = keyname;
				cfg.DataPtr = &(m->roi[i].mdl->par[IDX_TAU0 + 2 * j]);
				cfg.VarType = Cfg_Double;
				res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
				m->roi[i].mdl->par[IDX_TAU0 + 2 * j] *= 1.0 / m->tcal;
#ifdef DEBUG
				printf ("%s: %f\n", keyname,
						m->roi[i].mdl->par[IDX_TAU0 + 2 * j]);
#endif
				if (res < 1)
				{
					cfg.DataPtr = keyvalue;
					cfg.VarType = Cfg_String;
					ReadCfg (ini_file, sectionname, &cfg, CfgFile);
					set_link_fixed (m, IDX_TAU0 + 2 * j, keyvalue, i, 1);
				}
				else
				{
					set_variable (m, IDX_TAU0 + 2 * j, keyvalue, i);
				}

				sprintf (keyname, "amp[%d]", j);
				cfg.Name = keyname;
				cfg.DataPtr = &(m->roi[i].mdl->par[IDX_TAU0 + 2 * j + 1]);
				cfg.VarType = Cfg_Double;
				res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
#ifdef DEBUG
				printf ("%s: %f\n", keyname,
						m->roi[i].mdl->par[IDX_TAU0 + 2 * j + 1]);
#endif
				if (res < 1)
				{
					cfg.DataPtr = keyvalue;
					cfg.VarType = Cfg_String;
					ReadCfg (ini_file, sectionname, &cfg, CfgFile);
					set_link_fixed (m, IDX_TAU0 + 2 * j + 1, keyvalue, i, 0);
				}
				else if (((j + 1) == m->global_mdl->nexp)
						&& (m->global_mdl->fitType == LM_CHISQ_MLE))
				{
#ifdef DEBUG
					printf
					("last amplitude par[%d] given through normalization\n",
							IDX_TAU0 + 2 * j + 1);
#endif
				}
				else
				{
					set_variable (m, IDX_TAU0 + 2 * j + 1, keyvalue, i);
				}
			}
			// IRFSHIFT + BACKGROUND

			strcpy (keyname, "irfshift");
			cfg.Name = keyname;
			cfg.DataPtr = &(m->roi[i].mdl->par[IDX_IRFSHIFT]);
			cfg.VarType = Cfg_Double;
			res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
#ifdef DEBUG
			printf ("%s: %f\n", keyname, m->roi[i].mdl->par[IDX_IRFSHIFT]);
#endif
			if (res < 1)
			{
				cfg.DataPtr = keyvalue;
				cfg.VarType = Cfg_String;
				ReadCfg (ini_file, sectionname, &cfg, CfgFile);
				set_link_fixed (m, IDX_IRFSHIFT, keyvalue, i, 0);
			}
			else
			{
				set_variable (m, IDX_IRFSHIFT, keyvalue, i);
			}

			strcpy (keyname, "background");
			cfg.Name = keyname;
			cfg.DataPtr = &(m->roi[i].mdl->par[IDX_BG]);
			cfg.VarType = Cfg_Double;
			res = ReadCfg (ini_file, sectionname, &cfg, CfgFile);
#ifdef DEBUG
			printf ("%s: %f\n", keyname, m->roi[i].mdl->par[IDX_BG]);
#endif
			if (res < 1)
			{
				cfg.DataPtr = keyvalue;
				cfg.VarType = Cfg_String;
				ReadCfg (ini_file, sectionname, &cfg, CfgFile);
				set_link_fixed (m, IDX_BG, keyvalue, i, 0);
			}
			else
			{
				set_variable (m, IDX_BG, keyvalue, i);
			}
		}
	}

	create_pixellist (m);

	fclose (CfgFile);

	if (m->fileformat == PT3DATA) {
		if (!readpt3 (m->meas_fname, m->irf_fname, 0, m)){
			return NULL;
		}
	} else if (m->fileformat == QADATA) {
	    if (!readQA(m->meas_fname, m->irf_fname, 0, m)){
			return NULL;
		}
	}

	return m;
}



void
create_pixellist (FLIM_measurement * m)
{
	uint64_t i, j, checksum = 0;

	m->pixellist = calloc (m->width, sizeof (GSList **));
	m->pixelindex = calloc (m->width, sizeof (GSList **));

	for (i = 0; i < m->width; i++)
	{
		m->pixellist[i] = calloc (m->height, sizeof (GSList *));
		m->pixelindex[i] = calloc (m->height, sizeof (GSList *));
	}
	for (i = 0; i < m->nrois; i++)
	{
		for (j = 0; j < m->roi[i].npixels; j++)
		{
			// We have to save the corresponding ROI(s) [can be more than one] for each pixel
			m->pixellist[m->roi[i].pixel[j].x][m->roi[i].pixel[j].y] =
					g_slist_prepend (m->pixellist[m->roi[i].pixel[j].x]
					                              [m->roi[i].pixel[j].y], GUINT_TO_POINTER (i));


			// Wozu??
			m->pixelindex[m->roi[i].pixel[j].x][m->roi[i].pixel[j].y] =
					g_slist_prepend (m->pixelindex[m->roi[i].pixel[j].x]
					                               [m->roi[i].pixel[j].y], GUINT_TO_POINTER (j));
			checksum += m->roi[i].pixel[j].y + m->roi[i].pixel[j].x;
		}
#ifdef DEBUG
		printf ("checksum roi %ld: %ld\n", i, checksum);
#endif
		checksum = 0;
	}

	if (m->fileformat == QADATA) {
		m->irf_pixellist = calloc (m->width, sizeof (GSList **));
		m->irf_pixelindex = calloc (m->width, sizeof (GSList **));

		for (i = 0; i < m->width; i++)
		{
			m->irf_pixellist[i] = calloc (m->height, sizeof (GSList *));
			m->irf_pixelindex[i] = calloc (m->height, sizeof (GSList *));
		}
		for (i = 0; i < m->nirf_rois; i++)
		{
			for (j = 0; j < m->irf_roi[i].npixels; j++)
			{
				// We have to save the corresponding ROI(s) [can be more than one] for each pixel
				m->irf_pixellist[m->irf_roi[i].pixel[j].x][m->irf_roi[i].pixel[j].y] =
						g_slist_prepend (m->irf_pixellist[m->irf_roi[i].pixel[j].x]
						                              [m->irf_roi[i].pixel[j].y], GUINT_TO_POINTER (i));


				// Wozu??
				m->irf_pixelindex[m->irf_roi[i].pixel[j].x][m->irf_roi[i].pixel[j].y] =
						g_slist_prepend (m->irf_pixelindex[m->irf_roi[i].pixel[j].x]
						                               [m->irf_roi[i].pixel[j].y], GUINT_TO_POINTER (j));
				checksum += m->irf_roi[i].pixel[j].y + m->irf_roi[i].pixel[j].x;
			}
	#ifdef DEBUG
			printf ("checksum roi %ld: %ld\n", i, checksum);
	#endif
			checksum = 0;
		}
	}

	return;
}



int
calculate_roi_pixels (uint16_t xbottom, uint16_t ybottom, uint16_t xtop, uint16_t ytop,
		ROI *roi, FLIM_measurement *m)
{
	uint64_t j, k, l;

	if ((xbottom > m->width) || (ybottom > m->height))
	{
		printf ("Error in ROI definition:\n"
				"xbottom > width OR ybottom > height\n");
		printf ("width: %d, xbottom: %d\theight: %d, ybottom: %d\n", m->width, xbottom, m->height, ybottom);
		return 0;
	}

	if ((xtop >= xbottom) || (ytop >= ybottom))
	{
		printf ("Error in ROI definition:\n"
				"xtop >= xbottom OR ytop >= ybottom\n");
		return 0;
	}

	roi->npixels = (xbottom - xtop + 1) * (ybottom - ytop + 1);
#ifdef DEBUG
	printf ("npixels %ld\n", roi->npixels);
#endif
	roi->pixel = calloc (roi->npixels, sizeof (Pixel));
	l = 0;
	for (j = 0; j < (ybottom - ytop + 1); j++)
	{
		for (k = 0; k < (xbottom - xtop + 1); k++)
		{
			roi->pixel[l].photons = NULL;
			roi->pixel[l].x = xtop + k;
			roi->pixel[l].y = ytop + j;
			roi->pixel[l].ncounts = 0;
			roi->ncounts = 0;
			l++;
			//                              printf ("%d,%d\n", xtop + k,ytop + j);
		}
	}

	return 1;
}

int
calculate_roi_pixels_from_tiff_file (char *filename, FLIM_measurement *m) {
	TIFF *tif = TIFFOpen (filename, "r");
	int w,h,_h,x,y, roi, *pixel;
	uint32* raster;
	if (tif == NULL) {
		printf ("Could not open %s in order to read ROIs.\n", filename);
		return 0;
	}

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	if ((w != m->width) || (h != m->height)) {
		printf ("Dimensions of ROI containing file %s do not match the size specified in configuration file.\n", filename);
		return 0;
	}
	raster = _TIFFmalloc(w*h * sizeof (uint32));

	_h = h-1;

	if (raster != NULL) {
		if (TIFFReadRGBAImage(tif, w, h, raster, 0)) {
			m->nrois = 0;
			/* we need to loop several times over raster...
		this would be more compact if we would use linked lists for
		the rois...
			 */
			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
					if(((uint8) raster[(_h-y)*w+x]) != 0) {
						roi = (uint8) raster[(_h-y)*w+x];
						if (roi > m->nrois)
							m->nrois = roi;
					}
				}
			}
			if (m->nrois == 0) {
				printf ("Error while reading ROIs from %s: No rois found.\n", filename);
				return 0;
			}
			m->roi = create_rois (m->nrois, m->nsamples, m->global_mdl->nexp);
			m->joint_histogram = calloc (m->nrois * m->nsamples, sizeof (double));
#ifdef DEBUG
			printf ("nrois %d\n", m->nrois);
#endif

			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
					if(((uint8) raster[(_h-y)*w+x]) != 0) {
						roi = (uint8) raster[(_h-y)*w+x] - 1;
						m->roi[roi].npixels++;
					}
				}
			}
			for (roi = 0; roi < m->nrois; roi++) {
				m->roi[roi].pixel = calloc(m->roi[roi].npixels, sizeof (Pixel));
				m->roi[roi].ncounts = 0;
			}
			pixel = calloc(m->nrois,sizeof(uint32));
			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
					if(((uint8) raster[(_h-y)*w+x]) != 0) {
						roi = (uint8) raster[(_h-y)*w+x] - 1;
						m->roi[roi].pixel[pixel[roi]].x = x;
						m->roi[roi].pixel[pixel[roi]].y = y;
						m->roi[roi].pixel[pixel[roi]].ncounts = 0;
						m->roi[roi].pixel[pixel[roi]].photons = NULL;
						pixel[roi]++;
					}
				}
			}
			free(pixel);
		}
	}
	_TIFFfree(raster);
	TIFFClose(tif);
	return 1;
}

void
set_link_fixed (FLIM_measurement * m, uint16_t ipar, char *keyvalue_raw,
		uint16_t iroi, uint8_t islifetime)
{
	char *keyvalue;
	GList *ref_list;

	keyvalue = strtok (keyvalue_raw, "[]()");
	if (strcmp (keyvalue, "link") == 0)
	{
		uint16_t ilink, link_idx;
		Reference *ref, *hit_ref;
		GList *hit;
		GList *pvar_list;
		GList *link_pvar_list;
		uint16_t increment = 1;

		keyvalue = strtok (NULL, "[]");
		ilink = atoi (keyvalue);

		ref = malloc (sizeof (Reference));

		ref->iroi = iroi;
		ref->ipar = ipar;
		ref->ilink = ilink + 1;
		link_pvar_list = NULL;
		for (pvar_list = g_list_first(m->global_mdl->pvar_multi_map); pvar_list != NULL; pvar_list = pvar_list->next)
		{
			hit = g_list_first (pvar_list->data);	//, ref, (GCompareFunc) link_exists
			hit_ref = (Reference *) hit->data;
			if (hit_ref->ilink == ref->ilink)
			{
				increment = 0;
				link_pvar_list = pvar_list;
				break;
			}
		}

		if (increment)
		{
			m->global_mdl->npvar++;
			if (islifetime)
				m->global_mdl->link[ilink] *= 1.0 / m->tcal;
		}
		m->roi[iroi].mdl->par[ipar] = m->global_mdl->link[ilink];
		if (link_pvar_list != NULL) {
			ref_list = (GList *) link_pvar_list->data;
 			ref_list = g_list_append (ref_list, ref);
		}
		else {
			ref_list = NULL;
			ref_list = g_list_append (ref_list, ref);
			m->global_mdl->pvar_multi_map = g_list_append(m->global_mdl->pvar_multi_map,ref_list);
		}
		link_idx = g_list_index(m->global_mdl->pvar_multi_map, ref_list);
#ifdef DEBUG
		printf ("link_idx: %d\n",link_idx);
#endif
		m->global_mdl->pvar_multi[link_idx] = m->global_mdl->link[ilink];
#ifdef DEBUG
		printf ("par[%d] corresponds to link[%d]\n", ipar, atoi (keyvalue));
#endif
	}
	else if (strcmp (keyvalue, "fix") == 0)
	{
		keyvalue = strtok (NULL, "()");
		m->roi[iroi].mdl->par[ipar] = atof (keyvalue);
		if (islifetime)
			m->roi[iroi].mdl->par[ipar] *= 1.0 / m->tcal;
#ifdef DEBUG
		printf ("par[%d] fixed to %f\n", ipar, m->roi[iroi].mdl->par[ipar]);
#endif
	}

#ifdef DEBUG
	printf ("m->global_mdl->npvar: %d\n", m->global_mdl->npvar);
#endif
	return;
}


void
set_variable (FLIM_measurement * m, uint16_t ipar, char *keyvalue_raw, uint16_t iroi)
{
	Reference *ref;
	GList *ref_list;

	ref = malloc (sizeof (Reference));
	ref->iroi = iroi;
	ref->ipar = ipar;
	ref->ilink = 0;
	m->global_mdl->pvar_multi[m->global_mdl->npvar] =
			m->roi[iroi].mdl->par[ipar];

	ref_list = NULL;
	ref_list = g_list_append (ref_list, ref);
	m->global_mdl->pvar_multi_map = g_list_append(m->global_mdl->pvar_multi_map, ref_list);

	m->roi[iroi].mdl->pvar[m->roi[iroi].mdl->npvar] = ipar;
	m->roi[iroi].mdl->npvar++;
	m->global_mdl->npvar++;

#ifdef DEBUG
	printf ("m->global_mdl->npvar: %d\n", m->global_mdl->npvar);
#endif

	return;
}


void
read_globals_file (Histogram * measure, Model * conv_mdl, char *fname)
{
	uint16_t nsamples = 4096;
	double tcal;

	conv_mdl->nsamples = measure->nsamples = nsamples;
	measure->sample = malloc (sizeof (double) * nsamples);
	measure->irf = malloc (sizeof (double) * nsamples);



	if (checkfile (fname, &(conv_mdl->tstart), &(conv_mdl->tstop), &tcal) == 1)
	{
		measure->tcal = tcal;

		printf
		("\nLese Datei von Zeile %d bis Zeile %d, Zeitintervall: %.4f ns ein...",
				conv_mdl->tstart, conv_mdl->tstop, tcal);

		readfile (fname, conv_mdl->tstart,
				conv_mdl->tstop - conv_mdl->tstart, measure->sample,
				measure->irf);
		conv_mdl->irf_tstart = conv_mdl->tstart;
		conv_mdl->irf_tstop = conv_mdl->tstop;
		measure->sample += conv_mdl->tstart;
		measure->irf += conv_mdl->tstart;
		conv_mdl->nsamples = conv_mdl->tstop - conv_mdl->tstart;
		measure->nsamples = conv_mdl->tstop - conv_mdl->tstart;
		conv_mdl->irf = measure->irf;
		printf ("  done.\n");
	}

	return;
}



int
checkfile (char file[], uint16_t * rowstart, uint16_t * rowend, double *timestep)
{				//Identifiziere Start- und End-Zeile für den Import, ebenso die Größe des Zeitintervalls
	FILE *stream;
	char line[100];
	int firstline = 0;
	int success = 0;
	if ((stream = fopen (file, "r")) != NULL)
	{
		while (!feof (stream))
		{
			if (fgets (line, 100, stream) == NULL)
			{
				printf ("Datei wurde eingelesen.\n");
			}
			else
			{
				if (firstline == 0)
				{		//in der ersten Zeile: Zeitschritt, die Auswertung erfolgt zwischen den Punkten min und max
					sscanf (line, "%lf %hud %hu", timestep, rowstart, rowend);
					firstline = 1;
					success = 1;
					break;
				}
			}

		}
		fclose (stream);
	}
	else
	{
		printf ("Kann die Datei nicht oeffnen\n");
		success = 0;
	}
	return success;
}

void
readfile (char file[], uint64_t rowstart, uint64_t length,
		double *fluorescence, double *irf)
{
	FILE *stream;
	char line[100];
	int counter = 0;		//zeilenzähler

	if ((stream = fopen (file, "r")) != NULL)
	{
		while (!feof (stream))
		{
			if (fgets (line, 100, stream) == NULL)
			{
				printf ("Datei wurde eingelesen.\n");
			}
			else
			{
				if (counter > 0)
				{		//in der ersten Zeile: Zeitschritt, die Auswertung erfolgt zwischen den Punkten min und max
					sscanf (line, "%lf %lf", &irf[counter - 1],
							&fluorescence[counter - 1]);
				}
			}
			counter++;
		}
		fclose (stream);
	}
	else
	{
		printf ("Error opening file");
	}
	return;
}



void
write2file (double *array, uint64_t size)
{
	FILE *fp;

	printf ("Beginning to write file...\n");
	fp = fopen ("./test.bin", "wb");
	fwrite (array, sizeof (double), size, fp);
	fclose (fp);
	printf ("Done writing file...\n");

	return;
}

void
read_from_file (double *array, uint64_t size)
{
	FILE *fp;

	printf ("Opening file...\n");
	fp = fopen ("./test.bin", "rb");
	printf ("%ld elements read.\n", fread (array, sizeof (double), size, fp));
	fclose (fp);
	printf ("Done reading file...\n");

	return;
}


void
read_histogram_from_file (double *histogram, uint64_t total_size,
		char *fname, uint64_t start, uint64_t end)
{
	FILE *fp;
	uint64_t i;

	printf ("Reading histogram from file...\n");
	fp = fopen (fname, "r");
	for (i = 0; i < 11; i++)
	{
		while (fgetc (fp) != '\n');
	}
	for (i = 0; i < total_size; i++)
	{
		fscanf (fp, "%lf\n", &histogram[i]);
	}
	fclose (fp);

	memset (histogram, 0.0, start * sizeof (double));
	memset (&(histogram[end]), 0.0, (total_size - end) * sizeof (double));

	return;
}


void
write2expfile (double *amp, double *lifetime, double *meanamp,
		uint64_t meanampsize, double *ncounts, uint64_t ndifcounts,
		uint64_t nmodels)
{
	FILE *fp;
	uint64_t i, j;
	double *init_amp = amp;

	printf ("Beginning to write statistical data into file...\n");
	fp = fopen ("./2expamp.bin", "w");

	// write header
	for (j = 0; j < nmodels; j++)
	{
		for (i = 0; i < ndifcounts; i++)
		{
			double averagetau;

			averagetau = meanamp[i * nmodels + j] * lifetime[0] +
					(1.0 - meanamp[i * nmodels + j]) * lifetime[1];
			fprintf (fp, "%f %f %e\n", ncounts[i], amp[0],
					averagetau * 10.9e-12 / 10e-9);
		}
		fprintf (fp, "\n");
		amp += 2;
	}
	fclose (fp);
	printf ("Done writing statistical data into file...\n");
	amp = init_amp;

	return;
}


void
write2expampfixfile (double *variance, double variance_boundary,
		double *ncounts, uint64_t ndifcounts, double *meanamp)
{
	FILE *fp;
	uint64_t i;

	fp = fopen ("./2expampfix_neyman.csv", "w");

	printf ("Beginning to write statistical data into file...\n");
	fprintf (fp, "Counts\tsigma(alpha_1)\tsigma_cr(alpha_1)\tmeanamp\n");
	for (i = 0; i < ndifcounts; i++)
	{
		fprintf (fp, "%f\t%f\t%f\t%f\n", ncounts[i], sqrt (variance[i]),
				sqrt (1.0 / ncounts[i] * variance_boundary), meanamp[i]);
	}
	printf ("Done writing statistical data into file...\n");

	fclose (fp);

}

void
write_flim_image (double ***matrix, FLIM_measurement * m)
{
	TIFF *tif;
	char fname[255];
	uint64_t comp, x, y;
	float *buf;

	buf = malloc (sizeof(float)*m->width);
	sprintf (fname, "%s_amp.tiff", m->meas_fname);

	/* remove file, if exists */
	if (access( fname, F_OK ) != -1) {
		remove(fname);
	}

	for (comp = 0; comp < m->global_mdl->nexp; comp++)
	{
		tif = TIFFOpen (fname, "a");
	    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, m->width);
	    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, m->height);
	    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
	    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		for (y = 0; y < m->height; y++)
		{
			for (x = 0; x < m->width; x++)
			{
				buf[x] = matrix[x][y][comp];
			}
			TIFFWriteScanline (tif, buf, y, 0);
		}

		TIFFClose(tif);
	}
	free(buf);


	return;
}


void
write_flim_image_average_tau (double **matrix, FLIM_measurement * m)
{
	TIFF *tif;
	char fname[255];
	uint64_t x, y;
	float *buf;

	sprintf (fname, "%s_avg_tau.tiff",m->meas_fname);

	buf = _TIFFmalloc (sizeof(float)*m->width);

	tif = TIFFOpen (fname, "w");
	TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, m->width);
	TIFFSetField (tif, TIFFTAG_IMAGELENGTH, m->height);
	TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	for (y = 0; y < m->height; y++)
	{
		for (x = 0; x < m->width; x++)
		{
			buf[x] = matrix[x][y];
		}
		TIFFWriteScanline (tif, buf, y, 0);
	}
	TIFFClose(tif);
	_TIFFfree(buf);


	return;
}

void
write_levmar_result (FLIM_measurement * m, double *info)
{
	FILE *fp;
	double res;
	uint64_t sum_tchannels, i, j;
	double *function, res_roi[m->nrois][m->nsamples], chisq;

	fp = fopen ("global_results.csv", "w");

	fprintf (fp, "ROI\tParameter\tFitted value\n");
	for (i = 0; i < m->nrois; i++)
	{
		fprintf (fp, "%ld\t\"BG\"\t%f\n", i, m->roi[i].mdl->par[0]);
		fprintf (fp, "%ld\t\"IRFSHIFT\"\t%f\n", i, m->roi[i].mdl->par[1]);
		for (j = 0; j < m->global_mdl->nexp; j++)
		{
			fprintf (fp, "%ld\t\"Amp%ld\"\t%f\n", i, j,
					m->roi[i].mdl->par[IDX_TAU0 + 2 * j + 1]);
			fprintf (fp, "%ld\t\"Life%ld\"\t%f\n", i, j,
					m->roi[i].mdl->par[IDX_TAU0 + 2 * j] * m->tcal);
		}

	}
	// TODO parameter in einzelne rois übertragen
	for (i = 0; i < m->global_mdl->npvar; i++)
	{
		printf ("par[mdl->pvar[%ld]]: %f\n", i, m->global_mdl->pvar_multi[i]);
	}

	function = calloc (m->nsamples * m->nrois, sizeof (double));

	/* Convolve again because levmar has free()'d mdl->function */
	multi_convolve (m->global_mdl->pvar_multi, function,
			m->global_mdl->npvar, m->nsamples * m->nrois, m);

	sum_tchannels = 0;
	for (i = 0; i < m->nrois; i++)
	{
		uint64_t tchannels = m->roi[i].mdl->tstop - m->roi[i].mdl->tstart;

		chisq = compute_chisq_measure (res_roi[i], &(m->joint_histogram[m->nsamples * i]),	// &(m->roi[i].histogram[0]), //m->roi[i].mdl->tstart
				&(function[m->nsamples * i]),	// + m->roi[i].mdl->tstart
				m->nsamples, m->global_mdl->fitType, m);	//tchannels

		chisq *= 1.0 / (tchannels - m->roi[i].mdl->npvar);

		printf ("Chi square for roi[%ld]: %f\n", i, chisq);
		fprintf (fp, "%ld\t\"Chi square\"\t%f\n", i, chisq);

		sum_tchannels += tchannels;
	}
	printf ("Chi square: %f\nIterations: %f\n",
			info[1] / (sum_tchannels - m->global_mdl->npvar), info[5]);
	fprintf (fp, "-\t\"Global Chi square\"\t%f\n",
			info[1] / (sum_tchannels - m->global_mdl->npvar));
	printf ("FitType was %s.\n", fitType_str (m->global_mdl->fitType));
	fprintf (fp, "FitType was %s.\n", fitType_str (m->global_mdl->fitType));


	//      mdl->tstop-mdl->tstart
	printf ("Exit code: %f (%s)\n", info[6], lm_error_str (info[6]));
	fclose (fp);

	fp = fopen ("residuals.csv", "w");
	for (i = 0; i < m->nsamples * m->nrois; i++)
	{
		if (m->global_mdl->fitType == LM_CHISQ_MLE)
		{
			res = sqrt (function[i] - m->joint_histogram[i] -
					m->joint_histogram[i] * log (function[i] /
							m->joint_histogram[i]));
		}
		else
		{
			res =
					(function[i] -
							m->joint_histogram[i]) / sqrt (m->joint_histogram[i]);
		}
		fprintf (fp, "%f %f %f %f\n", m->joint_histogram[i], function[i],
				m->joint_histogram_irf[i], res);
	}
	fclose (fp);

	free (function);
}
