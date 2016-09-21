
typedef struct _QAData
{
 uint64_t NP;
 Photon * p;

} QAData;




QAData* NewQAData(uint64_t nphot);
int readQA(char *filename, char *filename_irf, uint64_t offset,
	   FLIM_measurement *m);
int writeQAData(char* filename,uint64_t offset,
	   FLIM_measurement *m, uint8_t *h);
int QAHeader(uint8_t* p, char* filename);
int PrintPhot(QAData* q);
