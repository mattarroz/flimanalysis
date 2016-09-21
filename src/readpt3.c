#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>

#include "pt3io.h"
#include "fit.h"
#include "readpt3.h"
#include "utils.h"


#define NO_IMG_HDR

int
readpt3 (char *filename, char *filename_irf,  uint64_t offset,
	 FLIM_measurement * m)
{

  uint64_t i, j, k;
  uint32_t datarate;
  FILE *fpin;
  float readtime;
  struct timeval tvstart, tvstop;
  GSList *iroi, *ipixel;
  tFullHdr *pt3hdr;
  tRecord Record;
  size_t result;
  uint64_t dlen;
  u_int64_t ofltime,outofrange,outofroi,outoftrange;
  double syncperiod,truensync, pixeldwelltime_sync;
  double linestart, framestart;
 #ifdef DOUBLE_MARKER4
  double dframe;
 #endif
  u_int16_t column; // Calculated by: (numsync-linestart)/m->pixeldwelltime
  u_int16_t row, frame;
  char betweenframes, betweenlines;
  
  if ((fpin = fopen (filename, "rb")) == NULL)
    {
      printf ("\ncannot open input file %s\n", filename);
      return 0;
    }

  pt3hdr = (tFullHdr *) malloc (sizeof (tFullHdr));
  if (!readpt3hdr (fpin, pt3hdr))
    {
      printf ("Error while reading pt3 header.\n");
      return 0;
    }


  dlen = 0;
  ofltime = 0;
  truensync = 0.0;
  row = column = frame = 0;
  framestart = 0.0;
  betweenframes = betweenlines = 1;
  outofrange = outofroi = outoftrange = 0;
  m->ncounts = 0;
  syncperiod = (double)1E9 / pt3hdr->TTTRHdr->CntRate0; //in nanoseconds
  pixeldwelltime_sync = m->pixeldwelltime / (syncperiod/1000.0); //m->pixeldwelltime has unit us, therefore /1000.0
  if (pt3hdr->TTTRHdr->ImgHdrSize == 0) {
	  pt3hdr->LSMImgHdr->PixX = m->width;
	  pt3hdr->LSMImgHdr->PixY = m->height;
  }

  gettimeofday (&tvstart, NULL);
  for (i = 0; (i < pt3hdr->TTTRHdr->Records) && (truensync < m->tabs_stop); i++)
  {
	  result = fread( &Record, 1, sizeof(tRecord) ,fpin);
	  if (result!= sizeof(tRecord))
	  {
		  printf("\nUnexpected end of input file!\n");
		  break;
	  }

	  if(Record.bits.channel==0xF) //this means we have a special record
	  {
		  switch (Record.special.markers) {
#ifdef NO_IMG_HDR
		  case 2:
#else
		  case 1: //Marker 1 means start of a line
#endif
			  if (betweenframes) break;
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  linestart = truensync;
			  column = 0;
			  betweenlines = 0;
			  break;
#ifdef NO_IMG_HDR
		  case 1:
#else
		  case 2: //Marker 2 means end of a line
#endif
			  if (betweenframes || betweenlines) break;
			  row++;
			  betweenlines = 1;
#ifndef DOUBLE_MARKER4
			  if (row == pt3hdr->LSMImgHdr->PixY) {
				  betweenframes = 1;
				  row = column = 0;
			  }
#endif
			  break;
		  case 4: //Marker 4 means end of frame if (truensync-framestart>=dframe)
#ifdef DOUBLE_MARKER4
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  if (frame == 0) framestart = truensync;
			  if (truensync-framestart>=dframe) {
				  framestart = truensync;
				  betweenframes = 1;
				  row = column = 0;
			  } else { // start of new frame
				  framestart = truensync;
				  betweenframes = 0;
				  frame++;
			  }
#else
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  framestart = truensync;
			  betweenframes = 0;
			  frame++;
#endif
			  break;
		  case 0: //not a marker means overflow
			  ofltime += T3WRAPAROUND; // unwrap the time tag overflow
			  break;
		  }

		  continue;
	  }



	  if(
			  (Record.bits.channel==0) //Should never occur in T3 Mode
			  ||(Record.bits.channel>4) //Should not occur with current routers
	  )
	  {
		  printf("\nIllegal Channel: #%1ld %1u\n",dlen,Record.bits.channel);
	  }

	  truensync = ((double)ofltime+(double)Record.bits.numsync);

	  if (!betweenlines && !betweenframes) {
		  column = (u_int16_t) round((truensync-linestart)/pixeldwelltime_sync);

		  if ((column < m->width) && (row < m->height)) {
			  if (m->pixellist[column][row])
			  {

				  for (iroi = m->pixellist[column][row], ipixel = m->pixelindex[column][row];
						  iroi;
						  iroi = g_slist_next (iroi), ipixel = g_slist_next (ipixel))
				  {
					  j = GPOINTER_TO_UINT (iroi->data);// printf ("j: %d ", j);
					  k = GPOINTER_TO_UINT (ipixel->data);

					  if ((Record.bits.dtime >= m->roi[j].mdl->tstart) && (Record.bits.dtime <= m->roi[j].mdl->tstop)) {
						  m->roi[j].histogram[Record.bits.dtime]++;
						  // The return value is the new start of the list
						  m->roi[j].pixel[k].photons =
								  g_slist_prepend (m->roi[j].pixel[k].photons,
										  GUINT_TO_POINTER (Record.bits.dtime));

						  m->roi[j].pixel[k].ncounts++;
						  m->roi[j].ncounts++;
						  m->ncounts++;
					  } else outoftrange++;
				  }
			  } else outofroi++;
		  } else outofrange++;

	  }

	  dlen++;
  }
  for (i = 0; i < m->nrois; i++)
  {
	  m->roi[i].mdl->ncounts = m->roi[i].ncounts;
  }

  fclose (fpin);

  gettimeofday (&tvstop, NULL);
  readtime =
    (tvstop.tv_sec + 1e-6 * tvstop.tv_usec) - (tvstart.tv_sec +
					       1e-6 * tvstart.tv_usec);
  datarate = (uint64_t) round (dlen / readtime);
  printf ("\nDone reading %ld photons from file %s in %f seconds (%d photons/s or %d bytes/s)...\n", dlen, filename, readtime, datarate, datarate * 14);	// + 1e-6*(float)tvstop.tv_usec + 1e-6*(float)tvstart.tv_usec
  free_pt3header(pt3hdr);

  printf ("%ld photons were outside specified image dimensions, %ld photons were outside the given ROIs, %ld photons were outside trange.\n", outofrange, outofroi, outoftrange);

  // Now read the Instrument Response Function (IRF)
  if ((fpin = fopen (filename_irf, "rb")) == NULL)
    {
      printf ("\ncannot open input file %s\n", filename_irf);
      return 0;
    }

  pt3hdr = (tFullHdr *) malloc (sizeof (tFullHdr));
  if (!readpt3hdr (fpin, pt3hdr))
    {
      printf ("Error while reading pt3 header.\n");
      return 0;
    }

  m->ncounts_irf = 0;
  dlen = 0;
  syncperiod = (double)1E9 / pt3hdr->TTTRHdr->CntRate0; //in nanoseconds
  gettimeofday (&tvstart, NULL);
  for (i = 0; (i < pt3hdr->TTTRHdr->Records) && (truensync < m->tabs_stop); i++)
  {
	  result = fread( &Record, 1, sizeof(tRecord) ,fpin);
	  if (result!= sizeof(tRecord))
	  {
		  printf("\nUnexpected end of input file!");
		  break;
	  }

	  if(Record.bits.channel==0xF) //this means we have a special record
	  {

		  // For confocal FLIM, we have no spatial resolution for the IRF,
		  //  therefore we need no to look for line/frame markers, as done
		  //  for the file read just before
		  if (Record.special.markers == 0) { //not a marker means overflow
			  ofltime += T3WRAPAROUND; // unwrap the time tag overflow
		  }

		  continue;
	  }


	  if(
			  (Record.bits.channel==0) //Should never occur in T3 Mode
			  ||(Record.bits.channel>4) //Should not occur with current routers
	  )
	  {
		  printf("\nIllegal Channel: #%1ld %1u",dlen,Record.bits.channel);
	  }

	  truensync = ((double)ofltime+(double)Record.bits.numsync);

	  m->joint_histogram_irf[Record.bits.dtime]++;
	  m->ncounts_irf++;

	  dlen++;
  }
  fclose (fpin);

  gettimeofday (&tvstop, NULL);
  readtime =
    (tvstop.tv_sec + 1e-6 * tvstop.tv_usec) - (tvstart.tv_sec +
					       1e-6 * tvstart.tv_usec);
  datarate = (uint64_t) round (m->ncounts_irf / readtime);
  printf ("\nDone reading %ld photons from file %s in %f seconds (%d photons/s or %d bytes/s)...\n",  m->ncounts_irf, filename_irf, readtime, datarate, datarate * 14);
  free_pt3header(pt3hdr);

  return 1;
}

int
readpt3_tabs (char *filename, char *filename_irf, uint64_t offset,
	 FLIM_measurement * m)
{

  uint64_t i, j, k;
  uint32_t datarate;
  FILE *fpin;
  float readtime;
  struct timeval tvstart, tvstop;
  GSList *iroi, *ipixel;
  tFullHdr *pt3hdr;
  tRecord Record;
  size_t result;
  uint64_t dlen;
  u_int64_t ofltime,outofrange,outofroi,outoftrange;
  double syncperiod,truensync, pixeldwelltime_sync;
  double linestart, framestart;
 #ifdef DOUBLE_MARKER4
  double dframe;
 #endif
  u_int16_t column; // Calculated by: (numsync-linestart)/m->pixeldwelltime
  u_int16_t row, frame;
  char betweenframes, betweenlines;

  if ((fpin = fopen (filename, "rb")) == NULL)
    {
      printf ("\ncannot open input file %s\n", filename);
      return 0;
    }

  pt3hdr = (tFullHdr *) malloc (sizeof (tFullHdr));
  if (!readpt3hdr (fpin, pt3hdr))
    {
      printf ("Error while reading pt3 header.\n");
      return 0;
    }


  dlen = 0;
  ofltime = 0;
  truensync = 0.0;
  row = column = frame = 0;
  framestart = 0.0;
  betweenframes = betweenlines = 1;
  outofrange = outofroi = outoftrange = 0;
  m->ncounts = 0;
  syncperiod = (double)1E9 / pt3hdr->TTTRHdr->CntRate0; //in nanoseconds
  pixeldwelltime_sync = m->pixeldwelltime / (syncperiod/1000.0); //m->pixeldwelltime has unit us, therefore /1000.0
  if (pt3hdr->TTTRHdr->ImgHdrSize == 0) {
	  pt3hdr->LSMImgHdr->PixX = m->width;
	  pt3hdr->LSMImgHdr->PixY = m->height;
  }

  gettimeofday (&tvstart, NULL);
  for (i = 0; (i < pt3hdr->TTTRHdr->Records) && (truensync < m->tabs_stop); i++)
  {
	  result = fread( &Record, 1, sizeof(tRecord) ,fpin);
	  if (result!= sizeof(tRecord))
	  {
		  printf("\nUnexpected end of input file!\n");
		  break;
	  }

	  if(Record.bits.channel==0xF) //this means we have a special record
	  {
		  switch (Record.special.markers) {
#ifdef NO_IMG_HDR
		  case 2:
#else
		  case 1: //Marker 1 means start of a line
#endif
			  if (betweenframes) break;
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  linestart = truensync;
			  column = 0;
			  betweenlines = 0;
			  break;
#ifdef NO_IMG_HDR
		  case 1:
#else
		  case 2: //Marker 2 means end of a line
#endif
			  if (betweenframes || betweenlines) break;
			  row++;
			  betweenlines = 1;
#ifndef DOUBLE_MARKER4
			  if (row == pt3hdr->LSMImgHdr->PixY) {
				  betweenframes = 1;
				  row = column = 0;
			  }
#endif
			  break;
		  case 4: //Marker 4 means end of frame if (truensync-framestart>=dframe)
#ifdef DOUBLE_MARKER4
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  if (frame == 0) framestart = truensync;
			  if (truensync-framestart>=dframe) {
				  framestart = truensync;
				  betweenframes = 1;
				  row = column = 0;
			  } else { // start of new frame
				  framestart = truensync;
				  betweenframes = 0;
				  frame++;
			  }
#else
			  truensync = ((double)ofltime+(double)Record.bits.numsync);
			  framestart = truensync;
			  betweenframes = 0;
			  frame++;
#endif
			  break;
		  case 0: //not a marker means overflow
			  ofltime += T3WRAPAROUND; // unwrap the time tag overflow
			  break;
		  }

		  continue;
	  }



	  if(
			  (Record.bits.channel==0) //Should never occur in T3 Mode
			  ||(Record.bits.channel>4) //Should not occur with current routers
	  )
	  {
		  printf("\nIllegal Channel: #%1ld %1u\n",dlen,Record.bits.channel);
	  }

	  truensync = ((double)ofltime+(double)Record.bits.numsync);

	  if (!betweenlines && !betweenframes) {
		  column = (u_int16_t) round((truensync-linestart)/pixeldwelltime_sync);

		  if (column < m->width) {
			  if (m->pixellist[column][row])
			  {

				  for (iroi = m->pixellist[column][row], ipixel = m->pixelindex[column][row];
						  iroi;
						  iroi = g_slist_next (iroi), ipixel = g_slist_next (ipixel))
				  {
					  j = GPOINTER_TO_UINT (iroi->data);// printf ("j: %d ", j);
					  k = GPOINTER_TO_UINT (ipixel->data);

					  if ((Record.bits.dtime >= m->roi[j].mdl->tstart) && (Record.bits.dtime <= m->roi[j].mdl->tstop)) {
						  m->roi[j].histogram[Record.bits.dtime]++;
						  // The return value is the new start of the list
						  m->roi[j].pixel[k].photons =
								  g_slist_prepend (m->roi[j].pixel[k].photons,
										  GUINT_TO_POINTER (Record.bits.dtime));

						  m->roi[j].pixel[k].ncounts++;
						  m->roi[j].ncounts++;
						  m->ncounts++;
					  } else outoftrange++;
				  }
			  } else outofroi++;
		  } else outofrange++;

	  }

	  dlen++;
  }
  for (i = 0; i < m->nrois; i++)
  {
	  m->roi[i].mdl->ncounts = m->roi[i].ncounts;
  }

  fclose (fpin);

  gettimeofday (&tvstop, NULL);
  readtime =
    (tvstop.tv_sec + 1e-6 * tvstop.tv_usec) - (tvstart.tv_sec +
					       1e-6 * tvstart.tv_usec);
  datarate = (uint64_t) round (dlen / readtime);
  printf ("\nDone reading %ld photons from file %s in %f seconds (%d photons/s or %d bytes/s)...\n", dlen, filename, readtime, datarate, datarate * 14);	// + 1e-6*(float)tvstop.tv_usec + 1e-6*(float)tvstart.tv_usec
  free_pt3header(pt3hdr);

  printf ("%ld photons were outside specified image dimensions, %ld photons were outside the given ROIs, %ld photons were outside trange.\n", outofrange, outofroi, outoftrange);

  // Now read the Instrument Response Function (IRF)
  if ((fpin = fopen (filename_irf, "rb")) == NULL)
    {
      printf ("\ncannot open input file %s\n", filename_irf);
      return 0;
    }

  pt3hdr = (tFullHdr *) malloc (sizeof (tFullHdr));
  if (!readpt3hdr (fpin, pt3hdr))
    {
      printf ("Error while reading pt3 header.\n");
      return 0;
    }

  m->ncounts_irf = 0;
  dlen = 0;
  syncperiod = (double)1E9 / pt3hdr->TTTRHdr->CntRate0; //in nanoseconds
  gettimeofday (&tvstart, NULL);
  for (i = 0; (i < pt3hdr->TTTRHdr->Records) && (truensync < m->tabs_stop); i++)
  {
	  result = fread( &Record, 1, sizeof(tRecord) ,fpin);
	  if (result!= sizeof(tRecord))
	  {
		  printf("\nUnexpected end of input file!");
		  break;
	  }

	  if(Record.bits.channel==0xF) //this means we have a special record
	  {

		  // For confocal FLIM, we have no spatial resolution for the IRF,
		  //  therefore we need no to look for line/frame markers, as done
		  //  for the file read just before
		  if (Record.special.markers == 0) { //not a marker means overflow
			  ofltime += T3WRAPAROUND; // unwrap the time tag overflow
		  }

		  continue;
	  }


	  if(
			  (Record.bits.channel==0) //Should never occur in T3 Mode
			  ||(Record.bits.channel>4) //Should not occur with current routers
	  )
	  {
		  printf("\nIllegal Channel: #%1ld %1u",dlen,Record.bits.channel);
	  }

	  truensync = ((double)ofltime+(double)Record.bits.numsync);

	  m->joint_histogram_irf[Record.bits.dtime]++;
	  m->ncounts_irf++;

	  dlen++;
  }
  fclose (fpin);

  gettimeofday (&tvstop, NULL);
  readtime =
    (tvstop.tv_sec + 1e-6 * tvstop.tv_usec) - (tvstart.tv_sec +
					       1e-6 * tvstart.tv_usec);
  datarate = (uint64_t) round (m->ncounts_irf / readtime);
  printf ("Done reading %ld photons from file %s in %f seconds (%d photons/s or %d bytes/s)...\n", m->ncounts_irf, filename_irf, readtime, datarate, datarate * 14);
  free_pt3header(pt3hdr);

  return 1;
}

