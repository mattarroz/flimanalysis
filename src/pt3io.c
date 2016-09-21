#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#define u_int64_t uint64_t
#define u_int32_t uint32_t
#define u_int16_t uint16_t
#define u_int8_t uint8_t
#include "pt3io.h"


tTxtHdr *read_pt3TxtHdr(FILE *fpin) {
 size_t result;
 tTxtHdr *TxtHdr = (tTxtHdr *) malloc(sizeof(tTxtHdr));

 result = fread( TxtHdr, 1, sizeof(tTxtHdr) ,fpin);
 if (result!= sizeof(tTxtHdr))
 {
  printf("\nerror reading txt header, aborted.");
  return NULL;
 }
 if( strcmp(TxtHdr->Ident,"PicoHarp 300") )
 {
    printf("\nFile identifier not found, aborted.");
	return NULL;
 }

 return TxtHdr;
}

int write_pt3TxtHdr(tTxtHdr *TxtHdr, FILE *fpout) {
 size_t result;

 result = fwrite( TxtHdr, 1, sizeof(tTxtHdr) ,fpout);
 if (result!= sizeof(tTxtHdr))
 {
  printf("\nerror writing txt header, aborted.");
  return 0;
 }
 return 1;
}

tTTTRHdr *read_pt3TTTRHdr(FILE *fpin) {
 size_t result;
 tTTTRHdr *TTTRHdr = (tTTTRHdr *) malloc(sizeof(tTTTRHdr));

 result = fread( TTTRHdr, 1, sizeof(tTTTRHdr) ,fpin);
 if (result!= sizeof(tTTTRHdr))
 {
  printf("\nerror reading TTTR header, aborted.");
  return NULL;
 }
	return TTTRHdr;
}

int write_pt3TTTRHdr(tTTTRHdr *TTTRHdr, FILE *fpout) {
 size_t result;
 
 result = fwrite( TTTRHdr, 1, sizeof(tTTTRHdr) ,fpout);
 if (result!= sizeof(tTTTRHdr))
 {
  printf("\nerror writing TTTR header, aborted.");
  return 1;
 }
	return 0;
}


tBoardHdr *read_pt3BoardHdr(FILE *fpin) {
 size_t result;
 tBoardHdr *BoardHdr = (tBoardHdr *) malloc(sizeof(tBoardHdr));

 result = fread( BoardHdr, 1, sizeof(tBoardHdr) ,fpin);
 if (result!= sizeof(tBoardHdr))
 {
  printf("\nerror reading Board header, aborted.");
  return NULL;
 }
	return BoardHdr;
}

int write_pt3BoardHdr( tBoardHdr *BoardHdr, FILE *fpout) {
 size_t result;
 
 result = fwrite( BoardHdr, 1, sizeof(tBoardHdr) ,fpout);
 if (result!= sizeof(tBoardHdr))
 {
  printf("\nerror writing Board header, aborted.");
  return 0;
 }
	return 1;
}


tBinHdr *read_pt3BinHdr(FILE *fpin) {
 size_t result;
 tBinHdr *BinHdr = (tBinHdr *) malloc(sizeof(tBinHdr));

 result = fread( BinHdr, 1, sizeof(tBinHdr) ,fpin);
 if (result!= sizeof(tBinHdr))
 {
  printf("\nerror reading Bin header, aborted.");
  return NULL;
 }
 if(BinHdr->MeasMode != MEASMODE_T3)
 {
    printf("\nWrong measurement mode, aborted.");
    return NULL;
 }
	return BinHdr;
}

int write_pt3BinHdr(tBinHdr *BinHdr, FILE *fpout) {
 size_t result;
 
 result = fwrite( BinHdr, 1, sizeof(tBinHdr) ,fpout);
 if (result!= sizeof(tBinHdr))
 {
  printf("\nerror writing Bin header, aborted.");
  return 0;
 }
	return 1;
}


tLSMImgHdr *read_pt3ImgHdr(FILE *fpin, tTTTRHdr *TTTRHdr) {
 /* skip the imaging header (you may need to read it if you
    want to interpret an imaging file */
 tLSMImgHdr *LSMImgHdr = (tLSMImgHdr *) malloc(sizeof(tLSMImgHdr));
 size_t result;

 if (TTTRHdr->ImgHdrSize == 0) {
	 printf ("Warning: TTTR file contains no image header.\n"
			"You should set all needed parameters manually.\n");
	 LSMImgHdr->Dimensions = 0;
	 LSMImgHdr->Ident = 0;
	 LSMImgHdr->Frame = 0;
	 LSMImgHdr->LineStart = 0;
	 LSMImgHdr->LineStop = 0;
	 LSMImgHdr->Pattern = 0;
	 LSMImgHdr->prcSin = 0;
	 LSMImgHdr->prcStart = 0;
	 LSMImgHdr->prcStop = 0;
	 LSMImgHdr->PixX = 0;
	 LSMImgHdr->PixY = 0;
	 LSMImgHdr->Resolution = 0;
	 return LSMImgHdr;
 }

 result = fread( LSMImgHdr, 1, sizeof(tLSMImgHdr) ,fpin);
 if (result!= sizeof(tLSMImgHdr))
 {
  printf("\nerror reading Image header, aborted.");
  return NULL;
 }
 
 fseek(fpin,TTTRHdr->ImgHdrSize*3,SEEK_CUR);
 
 return LSMImgHdr;
}

int write_pt3ImgHdr(tLSMImgHdr *LSMImgHdr, tTTTRHdr *TTTRHdr, FILE *fpout) {
 size_t result;

 result = fwrite( LSMImgHdr, 1, sizeof(tLSMImgHdr) ,fpout);
 if (result!= sizeof(tLSMImgHdr))
 {
  printf("\nerror writing Image header, aborted.");
  return 0;
 }
 
 fseek(fpout,TTTRHdr->ImgHdrSize*3,SEEK_CUR);
 
 return 1;
}

int
readpt3hdr (FILE * fpin, tFullHdr * pt3hdr)
{
  if ((pt3hdr->TxtHdr = read_pt3TxtHdr (fpin)) == NULL)
    return 0;
  if ((pt3hdr->BinHdr = read_pt3BinHdr (fpin)) == NULL)
    return 0;
  if ((pt3hdr->BoardHdr = read_pt3BoardHdr (fpin)) == NULL)
    return 0;
  if ((pt3hdr->TTTRHdr = read_pt3TTTRHdr (fpin)) == NULL)
    return 0;
  if ((pt3hdr->LSMImgHdr = read_pt3ImgHdr (fpin, pt3hdr->TTTRHdr)) == NULL)
    return 0;

  print_pt3header (pt3hdr->TxtHdr, pt3hdr->BinHdr, pt3hdr->BoardHdr, pt3hdr->TTTRHdr, pt3hdr->LSMImgHdr);

  return 1;
}


void free_pt3header(tFullHdr * pt3hdr) {
 free(pt3hdr->TxtHdr);
 free(pt3hdr->TTTRHdr);
 free(pt3hdr->BinHdr);
 free(pt3hdr->BoardHdr);
 free(pt3hdr->LSMImgHdr);
 free(pt3hdr);
}

void print_pt3header(tTxtHdr *TxtHdr, tBinHdr *BinHdr, tBoardHdr *BoardHdr, tTTTRHdr *TTTRHdr, tLSMImgHdr *LSMImgHdr){

 printf("Ident            : %.*s\n",(int)sizeof(TxtHdr->Ident),TxtHdr->Ident);
 printf("Format Version   : %.*s\n",(int)sizeof(TxtHdr->FormatVersion),TxtHdr->FormatVersion);
 printf("Creator Name     : %.*s\n",(int)sizeof(TxtHdr->CreatorName),TxtHdr->CreatorName);
 printf("Creator Version  : %.*s\n",(int)sizeof(TxtHdr->CreatorVersion),TxtHdr->CreatorVersion);
 printf("Time of Creation : %.*s\n",(int)sizeof(TxtHdr->FileTime),TxtHdr->FileTime);
 printf("File Comment     : %.*s\n",(int)sizeof(TxtHdr->CommentField),TxtHdr->CommentField);


 printf("No of Curves     : %d\n",BinHdr->Curves);
 printf("Bits per Record  : %d\n",BinHdr->BitsPerRecord);
 printf("RoutingChannels  : %d\n",BinHdr->RoutingChannels);
 printf("No of Boards     : %d\n",BinHdr->NumberOfBoards);
 printf("Active Curve     : %d\n",BinHdr->ActiveCurve);
 printf("Measurement Mode : %d\n",BinHdr->MeasMode);
 printf("Sub-Mode         : %d\n",BinHdr->SubMode);
 printf("Range No         : %d\n",BinHdr->RangeNo);
 printf("Offset           : %d\n",BinHdr->Offset);
 printf("AcquisitionTime  : %d\n",BinHdr->Tacq);
 printf("Stop at          : %d\n",BinHdr->StopAt);
 printf("Stop on Ovfl.    : %d\n",BinHdr->StopOnOvfl);
 printf("Restart          : %d\n",BinHdr->Restart);
 printf("DispLinLog       : %d\n",BinHdr->DispLinLog);
 printf("DispTimeAxisFrom : %d\n",BinHdr->DispTimeFrom);
 printf("DispTimeAxisTo   : %d\n",BinHdr->DispTimeTo);
 printf("DispCountAxisFrom: %d\n",BinHdr->DispCountsFrom);
 printf("DispCountAxisTo  : %d\n",BinHdr->DispCountsTo);





 printf("---------------------\n");

 printf(" HardwareIdent   : %.*s\n",(int)sizeof(BoardHdr->HardwareIdent),BoardHdr->HardwareIdent);
 printf(" HardwareVersion : %.*s\n",(int)sizeof(BoardHdr->HardwareVersion),BoardHdr->HardwareVersion);
 printf(" HardwareSerial  : %d\n",BoardHdr->HardwareSerial);
 printf(" SyncDivider     : %d\n",BoardHdr->SyncDivider);
 printf(" CFDZeroCross0   : %d\n",BoardHdr->CFDZeroCross0);
 printf(" CFDLevel0       : %d\n",BoardHdr->CFDLevel0 );
 printf(" CFDZeroCross1   : %d\n",BoardHdr->CFDZeroCross1);
 printf(" CFDLevel1       : %d\n",BoardHdr->CFDLevel1);
 printf(" Resolution/ns   : %lf\n",BoardHdr->Resolution);

 if(BoardHdr->RouterModelCode>0) //otherwise this information is meaningless
 {
   printf(" RouterModelCode       : %d\n",BoardHdr->RouterModelCode);  
   printf(" RouterEnabled         : %d\n",BoardHdr->RouterEnabled);   

   printf(" RtChan1_InputType     : %d\n",BoardHdr->RtChan1_InputType);
   printf(" RtChan1_InputLevel    : %d\n",BoardHdr->RtChan1_InputLevel); 
   printf(" RtChan1_InputEdge     : %d\n",BoardHdr->RtChan1_InputEdge);
   printf(" RtChan1_CFDPresent    : %d\n",BoardHdr->RtChan1_CFDPresent); 
   printf(" RtChan1_CFDLevel      : %d\n",BoardHdr->RtChan1_CFDLevel);
   printf(" RtChan1_CFDZeroCross  : %d\n",BoardHdr->RtChan1_CFDZeroCross);

   printf(" RtChan2_InputType     : %d\n",BoardHdr->RtChan2_InputType);
   printf(" RtChan2_InputLevel    : %d\n",BoardHdr->RtChan2_InputLevel); 
   printf(" RtChan2_InputEdge     : %d\n",BoardHdr->RtChan2_InputEdge);
   printf(" RtChan2_CFDPresent    : %d\n",BoardHdr->RtChan2_CFDPresent); 
   printf(" RtChan2_CFDLevel      : %d\n",BoardHdr->RtChan2_CFDLevel);
   printf(" RtChan2_CFDZeroCross  : %d\n",BoardHdr->RtChan2_CFDZeroCross);

   printf(" RtChan3_InputType     : %d\n",BoardHdr->RtChan3_InputType);
   printf(" RtChan3_InputLevel    : %d\n",BoardHdr->RtChan3_InputLevel); 
   printf(" RtChan3_InputEdge     : %d\n",BoardHdr->RtChan3_InputEdge);
   printf(" RtChan3_CFDPresent    : %d\n",BoardHdr->RtChan3_CFDPresent); 
   printf(" RtChan3_CFDLevel      : %d\n",BoardHdr->RtChan3_CFDLevel);
   printf(" RtChan3_CFDZeroCross  : %d\n",BoardHdr->RtChan3_CFDZeroCross);

   printf(" RtChan4_InputType     : %d\n",BoardHdr->RtChan4_InputType);
   printf(" RtChan4_InputLevel    : %d\n",BoardHdr->RtChan4_InputLevel); 
   printf(" RtChan4_InputEdge     : %d\n",BoardHdr->RtChan4_InputEdge);
   printf(" RtChan4_CFDPresent    : %d\n",BoardHdr->RtChan4_CFDPresent); 
   printf(" RtChan4_CFDLevel      : %d\n",BoardHdr->RtChan4_CFDLevel);
   printf(" RtChan4_CFDZeroCross  : %d\n",BoardHdr->RtChan4_CFDZeroCross);
 }



 printf("---------------------\n");



 printf("ExtDevices      : %d\n",TTTRHdr->ExtDevices);
 printf("CntRate0        : %d\n",TTTRHdr->CntRate0);
 printf("CntRate1        : %d\n",TTTRHdr->CntRate1);
 printf("StopAfter       : %d\n",TTTRHdr->StopAfter);
 printf("StopReason      : %d\n",TTTRHdr->StopReason);
 printf("Records         : %d\n",TTTRHdr->Records);
 printf("ImgHdrSize      : %d\n",TTTRHdr->ImgHdrSize);

 printf("---------------------\n");
 printf("Dimensions	 : %d\n", LSMImgHdr->Dimensions);
 printf("Ident		 : %d\n", LSMImgHdr->Ident);
 printf("Frame		 : %d\n", LSMImgHdr->Frame);
 printf("LineStart       : %d\n", LSMImgHdr->LineStart);
 printf("LineStop	 : %d\n", LSMImgHdr->LineStop);
 printf("Pattern	 : %s\n", LSMImgHdr->Pattern ? "Bidirectional" : "Monodirectional");
 printf("prcSin		 : %d\n", LSMImgHdr->prcSin);
 printf("prcStart	 : %d\n", LSMImgHdr->prcStart);
 printf("prcStop	 : %d\n", LSMImgHdr->prcStop);
 printf("PixX		 : %d\n", LSMImgHdr->PixX);
 printf("PixY		 : %d\n", LSMImgHdr->PixY);
 printf("Resolution:	 : %f\n", LSMImgHdr->Resolution);
 
printf("\nsync rate = %d / sec",TTTRHdr->CntRate0);
 printf("\nsync period = %5.4lf ns",(double)1E9 / TTTRHdr->CntRate0);
 printf("\n%d records",TTTRHdr->Records);


 return;
}
