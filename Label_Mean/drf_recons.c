#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
// #include "psf.h"
// #include "global_var_dec.h"


float LabelMean(float *pImgPre, float *pImgMean, float *pImgLabel, int length,int NumLabel)
{
   // time_t starttime,stoptime;
	clock_t tin,tout;
	float *label=calloc(NumLabel, sizeof(float));
	float *mean=calloc(NumLabel, sizeof(float));
	int *count=calloc(NumLabel, sizeof(int));
	int i,j;
	//starttime=time(NULL);
        tin=clock();
        printf("Do LabelMean operation\n");
	for (i=0;i<NumLabel;i++)
	{	
		label[i]=0/*i+1*/;
		mean[i]=0;
		count[i]=0;
	}
	for (i=0;i<length;i++)
	{	
		for (j=0;j<NumLabel;j++)
			if (pImgLabel[i]==label[j])
			{	mean[j]=mean[j]+pImgPre[i];
				count[j]++;
			}
	}
	for (j=0;j<NumLabel;j++)
		{mean[j]=mean[j]/count[j];	printf(" %f",mean[j]);}

	for (i=0;i<length;i++)
		for (j=0;j<NumLabel;j++)
			if (pImgLabel[i]==label[j])
			{
				pImgMean[i]=mean[j];
				// printf(" %f",mean[j]);
				break;
			}
	free(label);
	free(mean);
	free(count);	
	//stoptime=time(NULL);
	tout=clock();
	return (float)(tout-tin)/(float)CLOCKS_PER_SEC;//stoptime-starttime;
}


int main()
{
	float *pImgPre = (float*) calloc(288*416*120,sizeof(float));
	float *pImgMean = (float*) calloc(288*416*120,sizeof(float));
	float *pImgLabel = (float*) calloc(288*416*120,sizeof(float));
    memset(pImgPre,0,288*416*120*sizeof(float));
    memset(pImgMean,0,288*416*120*sizeof(float));
    memset(pImgLabel,0,288*416*120*sizeof(float));

    char out_name[300];
	int NumLabel = 3;
	int length = 288*416*120;
	FILE *iImgPre;
    iImgPre = fopen("../Test2016April/GPU_DHAPET_backward_ite1_attenuation_r001_new.DAT","rb");
    fread(pImgPre,sizeof(float),length,iImgPre);
    fclose(iImgPre);

    FILE *iImgLabel;
    iImgLabel = fopen("ImgLabel","rb");
    fread(pImgLabel,sizeof(float),length,iImgLabel);
    fclose(iImgLabel);

    LabelMean(pImgPre, pImgMean, pImgLabel, length, NumLabel);

    printf("Save data to ImgMean\n");
    sprintf(out_name,"ImgMean2.DAT");
    FILE *iImgMean;
    if((iImgMean=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
        exit(0);
    }
    fwrite(pImgMean,sizeof(float),length,iImgMean);
    fclose(iImgMean);

    // sprintf(out_name,"ImgPre.DAT");
    // FILE *iImgPre_1;
    // if((iImgPre_1=fopen(out_name,"wb"))==NULL)
    // {
    //     printf("Can not open file\n");
    //     exit(0);
    // }
    // fwrite(pImgPre,sizeof(float),length,iImgPre_1);
    // fclose(iImgPre_1);


}