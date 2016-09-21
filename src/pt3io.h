
#define DISPCURVES 8
#define T3WRAPAROUND 65536
#define MEASMODE_T2 2 
#define MEASMODE_T3 3


/*
The following structures are used to hold the file data
They directly reflect the file structure.
The data types used here to match the file structure are correct
for the tested compilers.
They may have to be changed for other compilers.
*/


#pragma pack(4) //structure alignment to 4 byte boundaries

/* These are substructures used below */

typedef struct _tParamStruct { float Start;
                float Step;
				float End;  } tParamStruct;

typedef struct _tCurveMapping{ int32_t MapTo;
				int32_t Show; } tCurveMapping;

/* The following represents the readable ASCII file header portion */

typedef struct _tTxtHdr {		char Ident[16];				//"PicoHarp 300"
				char FormatVersion[6];		//file format version
				char CreatorName[18];		//name of creating software
				char CreatorVersion[12];	//version of creating software
				char FileTime[18];
				char CRLF[2];
				char CommentField[256]; } tTxtHdr;

/* The following is binary file header information */

typedef struct _tBinHdr {		int32_t Curves;
				int32_t BitsPerRecord;
				int32_t RoutingChannels;
				int32_t NumberOfBoards;
				int32_t ActiveCurve;
				int32_t MeasMode;
				int32_t SubMode;
				int32_t RangeNo;
				int32_t Offset;
				int32_t Tacq;				// in ms
				int32_t StopAt;
				int32_t StopOnOvfl;
				int32_t Restart;
				int32_t DispLinLog;
				int32_t DispTimeFrom;		// 1ns steps
				int32_t DispTimeTo;
				int32_t DispCountsFrom;
				int32_t DispCountsTo;
				tCurveMapping DispCurves[DISPCURVES];	
				tParamStruct Params[3];
				int32_t RepeatMode;
				int32_t RepeatsPerCurve;
				int32_t RepeatTime;
				int32_t RepeatWaitTime;
				char ScriptName[20];	} tBinHdr;

/* The next is a board specific header */

typedef struct _tBoardHdr {		
				char HardwareIdent[16]; 
				char HardwareVersion[8]; 
				int32_t HardwareSerial; 
				int32_t SyncDivider;
				int32_t CFDZeroCross0;
				int32_t CFDLevel0;
				int32_t CFDZeroCross1;
				int32_t CFDLevel1;
				float Resolution;
				//below is new in format version 2.0
				int32_t RouterModelCode;
				int32_t RouterEnabled;
				int32_t RtChan1_InputType; 
				int32_t RtChan1_InputLevel;
				int32_t RtChan1_InputEdge;
				int32_t RtChan1_CFDPresent;
				int32_t RtChan1_CFDLevel;
				int32_t RtChan1_CFDZeroCross;
				int32_t RtChan2_InputType; 
				int32_t RtChan2_InputLevel;
				int32_t RtChan2_InputEdge;
				int32_t RtChan2_CFDPresent;
				int32_t RtChan2_CFDLevel;
				int32_t RtChan2_CFDZeroCross;
				int32_t RtChan3_InputType; 
				int32_t RtChan3_InputLevel;
				int32_t RtChan3_InputEdge;
				int32_t RtChan3_CFDPresent;
				int32_t RtChan3_CFDLevel;
				int32_t RtChan3_CFDZeroCross;
				int32_t RtChan4_InputType; 
				int32_t RtChan4_InputLevel;
				int32_t RtChan4_InputEdge;
				int32_t RtChan4_CFDPresent;
				int32_t RtChan4_CFDLevel;
				int32_t RtChan4_CFDZeroCross;		} tBoardHdr;


/* The next is a TTTR mode specific header */

typedef struct _tTTTRHdr{
				int32_t ExtDevices;
				int32_t Reserved1;
				int32_t Reserved2;			
				int32_t CntRate0;
				int32_t CntRate1;
				int32_t StopAfter;
				int32_t StopReason;
				int32_t Records;
				int32_t ImgHdrSize;		} tTTTRHdr;

typedef struct _tE710ImgHdr { int32_t Dimensions;
				int32_t Ident;
				int32_t TimePerPixel;
				int32_t Acceleration;
				int32_t Pattern;
				int32_t Reserved;
				float X0;
				float Y0;
				int32_t PixX;
				int32_t PixY;
				float PixResol;
				float TStartTo;
				float TStopTo;
				float TStartFro;
				float TStopFro;
} tE710ImgHdr;





typedef struct _tKDT180ImgHdr { int32_t Dimensions;
				int32_t Ident;
				int32_t Velocity;
				int32_t Acceleration;
				int32_t Pattern;
				int32_t Reserved;
				float X0;
				float Y0;
				int32_t PixX;
				int32_t PixY;
				float PixResol;
} tKDT180ImgHdr;


typedef struct _tLSMImgHdr{	int32_t Dimensions; //0, 1=reserved, 2=line, 3=area
				int32_t Ident; //3= LSM
				int32_t Frame; //frame trigger index
				int32_t LineStart; //line start trigger index
				int32_t LineStop; //line stop trigger index
				int8_t Pattern; //0=monodirectional scan; 1=bidirectional scan
				int8_t prcSin; //(reserved for measurements via TCP/IP-protocol)
				int8_t prcStart; //(reserved for measurements via TCP/IP-protocol)
				int8_t prcStop; //(reserved for measurements via TCP/IP-protocol)
				int32_t PixX; //number of pixels (X)
				int32_t PixY; //number of pixels (Y)
				float Resolution; //(reserved for measurements via TCP/IP-protocol)
} tLSMImgHdr;



/*The following data records appear for each T3 mode event*/

typedef union _tRecord	{ 
				 u_int32_t allbits;
				 struct
				 {	
					unsigned numsync	:16; 
					unsigned dtime		:12; 		
					unsigned channel	:4;
				 } bits;
				 struct
				 {	
					unsigned numsync	:16; 
					unsigned markers	:12;  					
					unsigned channel	:4;
				 } special;
		} tRecord;


typedef struct _tFullHdr {
	tTxtHdr *TxtHdr;
	tBinHdr *BinHdr;
	tBoardHdr *BoardHdr;
	tTTTRHdr *TTTRHdr;
	tLSMImgHdr *LSMImgHdr;
} tFullHdr;

int
readpt3hdr (FILE * fpin, tFullHdr * pt3hdr);
tTxtHdr *read_pt3TxtHdr(FILE *fpin);
int write_pt3TxtHdr(tTxtHdr *TxtHdr, FILE *fpout);
tBinHdr *read_pt3BinHdr(FILE *fpin);
int write_pt3BinHdr(tBinHdr *BinHdr, FILE *fpout);
tBoardHdr *read_pt3BoardHdr(FILE *fpin);
int write_pt3BoardHdr( tBoardHdr *BoardHdr, FILE *fpout);
tTTTRHdr *read_pt3TTTRHdr(FILE *fpin);
int write_pt3TTTRHdr(tTTTRHdr *TTTRHdr, FILE *fpout);
tLSMImgHdr *read_pt3ImgHdr(FILE *fpin, tTTTRHdr *TTTRHdr);
int write_pt3ImgHdr(tLSMImgHdr *LSMImgHdr, tTTTRHdr *TTTRHdr, FILE *fpout);
void print_pt3header(tTxtHdr *TxtHdr, tBinHdr *BinHdr, tBoardHdr *BoardHdr, tTTTRHdr *TTTRHdr, tLSMImgHdr *LSMImgHdr);
void free_pt3header(tFullHdr *pt3hdr);
