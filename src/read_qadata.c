#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>


#include "fit.h"
#include "read_qadata.h"
#include "utils.h"



int readQA(char *filename, char *filename_irf, uint64_t offset,
		FLIM_measurement *m)
{

	uint64_t i, j, k, nphotons;
	uint16_t W = m->width, H = m->height;
	uint64_t datarate;
	u_int64_t outofrange,outofroi,outoftrange;
	FILE *f;
	uint8_t *p, h[256];
	//    UINT *pixelmap;
	float readtime;
	uint16_t x, y, TAC;//, *reTAC;//, *photons,
	double tabs = 0.0, tabsold = 0.0, overflow = 0.0;
	struct timeval tvstart, tvstop;
	struct stat st;
	GSList *iroi, *ipixel;

	if (access( filename, F_OK ) == -1) {
		printf("Error while reading QAData file %s\n", filename);
		return 0;
	}
	stat(filename, &st);
	nphotons = (st.st_size - 256)/14;
	m->ncounts = 0;

	gettimeofday (&tvstart, NULL);
	printf ("Reading from file %s (size: %ld photons) ...\n", filename, nphotons);


	f = fopen(filename, "r+b");
	if (NULL == f) {
		printf("File '%s' not found", filename);
		return 0;
	}
	if (!fread(h, 1, sizeof(uint8_t) * 256, f)) {
		printf("Error: Unexpected end of file.\n");
		return 0;
	}

	fseek(f, 256 + offset * 14, SEEK_SET);
	p = (uint8_t *) calloc(1, 14);

	tabs = 0.0;
	for (i = 0; (i < nphotons) && (tabs < m->tabs_start); i++) {
		if (!fread(p, 1, 14, f)) {
			printf("Error: Unexpected end of file.\n");
			return 0;
		}
		tabs = (*(uint32_t*) (p + 10)/m->tabs_clock);

		if (tabs-tabsold < -2000.0)
		{
			printf ("Overflow. TAbs: %f, last TAbs: %f\n", tabs,
					tabsold);
			overflow += tabsold;
		}
		tabsold = tabs;
		tabs += overflow;

	}

	outofrange = 0;
	outoftrange = 0;
	outofroi = 0; 
	for (i = 0; (i < nphotons) && (tabs < m->tabs_stop); i++) {
		if (!fread(p, 1, 14, f)) {
			printf("Error: Unexpected end of file.\n");
			return 0;
		}
		y = round(*(float *) p * (H - 1));
		x = round(*(float *) (p + 4) * (W - 1));
		tabs = *(uint32_t*) (p + 10)/m->tabs_clock;

		if (tabs-tabsold < -2000.0)
		{
			printf ("Overflow. TAbs: %f, last TAbs: %f\n", tabs,
					tabsold);
			overflow += tabsold;
		}
		tabsold = tabs;
		tabs += overflow;

		//	printf ("x %d\t, y %d\n", x, y);
		if ((x < m->width) && (y < m->height)) {
			if(m->pixellist[x][y]) { 
				/* TAC is the bin number in the lifetime histogram */
				/* We need to reverse the histogram for QA detectors,
				 * therefore we subtract it from m->nsamples
				 * - 1 appears because m->nsamples is a size
				 * 	but we want to use TAC as array index
				 */ 
				TAC = m->nsamples - *(uint16_t *) (p + 8) - 1;

				for (iroi = m->pixellist[x][y], ipixel = m->pixelindex[x][y];
						iroi;
						iroi = g_slist_next (iroi), ipixel = g_slist_next(ipixel)) {
					j = GPOINTER_TO_UINT (iroi->data);
					k = GPOINTER_TO_UINT (ipixel->data);
					if (!((TAC >= m->roi[j].mdl->tstart) && (TAC <= m->roi[j].mdl->tstop))) {
						outoftrange++;
						continue;
					}

					m->joint_histogram[j*m->nsamples + TAC]++;
					m->roi[j].pixel[k].photons =
							g_slist_prepend(m->roi[j].pixel[k].photons,
									GUINT_TO_POINTER(TAC));

					m->roi[j].pixel[k].ncounts++;
					m->roi[j].ncounts++;
					m->ncounts++;
				}
			} else outofroi++;
		} else outofrange++;
	}
	for (i = 0; i < m->nrois; i++) {
		m->roi[i].mdl->ncounts = m->roi[i].ncounts;
	} 

	fclose(f);


	gettimeofday (&tvstop, NULL);
	readtime = (tvstop.tv_sec + 1e-6*tvstop.tv_usec) - (tvstart.tv_sec + 1e-6*tvstart.tv_usec);
	datarate =  nphotons / readtime;
	
	printf ("Done reading %lld photons from file %s in %f seconds (%lld photons/s or %lld bytes/s)...\n",
			 nphotons,filename, readtime, datarate, datarate*14); // + 1e-6*(float)tvstop.tv_usec + 1e-6*(float)tvstart.tv_usec


	printf ("%ld photons were outside specified image dimensions, %ld photons were outside the given ROIs, %ld photons were outside trange.\n", outofrange, outofroi, outoftrange);

	//    writeQAData("qatest2.QAData", 0, m, h);


	if (access( filename_irf, F_OK ) == -1) {
		printf("Error while reading QAData file %s\n", filename_irf);
		return 0;
	}
	stat(filename_irf, &st);
	nphotons = (st.st_size - 256)/14;
	m->ncounts_irf = 0;

	f = fopen(filename_irf, "r+b");
	if (NULL == f) {
		printf("File '%s' not found", filename);
		return 0;
	}
	if (!fread(h, 1, sizeof(uint8_t) * 256, f)) {
		printf("Error: Unexpected end of file.\n");
		return 0;
	}

	fseek(f, 256 + offset * 14, SEEK_SET);

	gettimeofday (&tvstart, NULL);
	printf ("Reading %lld photons from file %s ...\n", nphotons, filename_irf);

	outofrange = 0;
	outoftrange = 0;
	outofroi = 0;
	for (i = 0; i < nphotons; i++) {
		if (!fread(p, 1, 14, f)) {
			printf("Error: Unexpected end of file.\n");
			return 0;
		}
		x = round(*(float *) p * (W - 1));
		y = round(*(float *) (p + 4) * (H - 1));

		//	printf ("x %d\t, y %d\n", x, y);
		if ((x < m->width) && (y < m->height)) {
			if(m->irf_pixellist[x][y]) {
				TAC = m->nsamples - *(uint16_t *) (p + 8) - 1;

				for (iroi = m->irf_pixellist[x][y], ipixel = m->irf_pixelindex[x][y];
						iroi;
						iroi = g_slist_next (iroi), ipixel = g_slist_next(ipixel)) {
					j = GPOINTER_TO_UINT (iroi->data);
					k = GPOINTER_TO_UINT (ipixel->data);
					if (!((TAC >= m->irf_roi[j].tstart) && (TAC <= m->irf_roi[j].tstop))) {
						outoftrange++;
						continue;
					}

					m->joint_histogram_irf[j*m->nsamples + TAC]++;
					m->irf_roi[j].pixel[k].ncounts++;
					m->irf_roi[j].ncounts++;
					m->ncounts_irf++;
				}
			} else outofroi++;
		} else outofrange++;
	}

	gettimeofday (&tvstop, NULL);
	readtime = (tvstop.tv_sec + 1e-6*tvstop.tv_usec) - (tvstart.tv_sec + 1e-6*tvstart.tv_usec);
	datarate =  nphotons / readtime;
	
	printf ("Done reading %lld photons from file %s in %f seconds (%lld photons/s or %lld bytes/s)...\n",
			 nphotons,filename_irf, readtime, datarate, datarate*14); // + 1e-6*(float)tvstop.tv_usec + 1e-6*(float)tvstart.tv_usec

	
	printf ("%ld photons were outside specified image dimensions, %ld photons were outside the given ROIs, %ld photons were outside trange.\n", outofrange, outofroi, outoftrange);


	free(p);
	return 1;
}

int writeQAData(char* filename,uint64_t offset,
		FLIM_measurement *m, uint8_t *h)
{
	FILE* f;
	uint64_t i, j, counter = 256;
	uint32_t TAC = 0;
	uint64_t TAbs = 0;
	GSList *photon;
	float x, y;

	f=fopen(filename,"wb");
	fwrite(h,256,1,f);
	for(i = 0; i < m->nrois; i++)
	{
		for (j = 0; j < m->roi[i].npixels; j++) {
			for (photon = m->roi[i].pixel[j].photons; //m->roi[i].pixel[j].first_photon
					photon; photon = g_slist_next(photon)) {
				x = (float)m->roi[i].pixel[j].x/(float)m->width;
				y = (float)m->roi[i].pixel[j].y/(float)m->height;
				fwrite(&x,1,4,f);
				//		     printf ("x %f, m->roi[i].pixel[j].x %d, m->width %d\n", x, m->roi[i].pixel[j].x, m->width);
				fwrite(&y,1,4,f);
				//		     printf ("y %f, m->roi[i].pixel[j].y %d, m->height %d\n", y, m->roi[i].pixel[j].y, m->height);
				TAC = GPOINTER_TO_UINT(photon->data);
				fwrite(&TAC,1,2,f);
				//		     printf ("TAC %d\n", TAC);
				TAbs = counter*1000;
				fwrite(&TAbs,1,4,f);
				//		     printf ("Tabs %ld\n", TAbs);
				counter += 14;
			}
		}
	}
	printf ("counter %ld\n", counter);

	fclose(f);
	return 1;
}


