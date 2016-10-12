#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#define _GNU_SOURCE
#include <getopt.h>
#define Tr = 100
#include "global_var_dec.h"
#include "psf_analytic.h"

float s_time()
{
 return (float)clock()/CLOCKS_PER_SEC;
}

int nboundary(int cy, int cz)
{
    if( cy < ncy_k && cy >= 0 &&
        cz < ncz_k && cz >= 0    )
        return 1;
    else return 0;
}

int boundary(int cy, int cz)
{
    if( cy < ncy && cy >= 0 &&
        cz < ncz && cz >= 0    )
        return 1;
    else return 0;
}

int boundary_z(int cz)
{
    if(cz < ncz && cz >= 0) return 1;
    else return 0;
}

int boundary_y(int cy)
{
    if(cy < ncy && cy >= 0) return 1;
    else return 0;
}

void config()
{
    nkernel=16;
    nimgy = 4 * ncy;
    nimgz = 4 * ncz;    
    ncrow = ncy * ncz;
    nccol = ncy * ncz;
    ncy_k = 2 * ncz;
    ncz_k = 2 * ncz;
    ncrow_k = ncy_k * ncz_k;
    nccol_k = ncrow_k;

    printf("  ncy=%d, ncz=%d\n",ncy,ncz); 
    printf("  pitch=%.2f\n",pitch); 
    printf("  nimgy=%d, nimgz=%d, nimgx=%d\n",nimgy,nimgz,nimgx); 
    printf("  ngrid=%d\n",ngrid); 
    printf("  spacing=%2.0fmm\n",spacing);
    printf("  index of xplane=%d\n",x_plane);
    //printf("  the object data will be saved with description: %s\n",description_name);
    py_edge_low = -ncy_k * cw / 2.;
    py_edge_up  =  ncy_k * cw / 2.;
    pz_edge_low = -ncz_k * cw / 2.;
    pz_edge_up  =  ncz_k * cw / 2.;
    vx = spacing / nimgx;
    vy = cw / (4);
    vz = cw / (4);

    kvx_edge_low = vx * (x_plane - 0.5);
    kvx_edge_up  = vx * (x_plane + 0.5);

}

void forward(float* pFwd,float* pBwd, int t, int tof)
{
 
 int index;
 float value;
 int count=0;
 int fx_plane = x_plane + 60;
 //x_plane=-1;
 TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));

 int i,j;

 FILE *inTAB,*inINCRE,*inSparse,*intSparse;


 for(i=0;i<4;i++)
 {
  for(j=0;j<4;j++)
  {
   char name[500];
   //char kernel_name[100][10];
   sprintf(kernel_name[i*4+j],"y%dz%d",j,i);

   sprintf(name, "../data_tof/sparse_table_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
   int* tSparse= (int*) calloc(pTable[nccol_k].bias,sizeof(int));  

   sprintf(name, "../data_tof/sparse_incre_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../data_tof/sparse_matrix_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   inSparse = fopen(name,"rb");
   fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse); 

   sprintf(name, "../data_tof/sparse_tof_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   intSparse = fopen(name,"rb");
   fread(tSparse,sizeof(int),pTable[nccol_k].bias,intSparse); 


   int LOR_initial=pTable[nccol_k].bias;

   int iz,jy,kcz,kcy;
   for(iz=0;iz<ncz;iz++)
   {
   
    int shift_z= (iz*4-ncz*2)/4;
    for(jy=0;jy<ncy;jy++)
    {
     
     int shift_y = (jy*4-ncy*2)/4;  
     for(kcz=0;kcz<ncz;kcz++)
     {
      int cz1=kcz+52-shift_z;
      for(kcy=0;kcy<ncy;kcy++)
      {
       int cy1=kcy+68-shift_y;
       if(nboundary(cy1,cz1))
       {
        int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
        int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
        int c2total=pTable[cz1*ncy_k+cy1+1].bias;
  	int ktotal;
        for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
        {
         int cy2=pIncre[ktotal]._y+cy2min;
         int cz2=pIncre[ktotal]._z+cz2min;
	 if(boundary(cy2-68,cz2-52))
         {
	 
	  if(tSparse[ktotal]==t)
	  {
          		pFwd[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal]*pBwd[fx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]/(2880*720);
			count+=1;
			

	  }
	  
         }
        }
       }
      }
     }
    }
   }
   fclose(inTAB);
   fclose(inINCRE);
   fclose(inSparse);
   fclose(intSparse);
   free(pIncre);
   free(pSparse);
   free(tSparse);
  }
 }
 free(pTable);
 printf("tof=%d,count=%d\n",tof,count);
} 

void ratio(float *pFW, float *pHist, float *pRatio, int ite, float *cost_function)
{
    int forcount=0;
    int hiscount=0;
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {
	if(pFW[i]>0)
	{
		forcount+=1;
	}
	if(pHist[i] >= 0)
	{
		hiscount+=1;
	}
        if(pHist[i] >= 0)
        {
	    
            if(pFW[i] > 0)
            {
                pRatio[i] =  pHist[i] / pFW[i];
		
            }
            else
            {
                if(pHist[i] > 0)
                {
                    pRatio[i] = 0; // non_zero / zero
                }
                else
                { 
                    pRatio[i] = 1.0; // zero / zero
                }
            }
        }
        if(pRatio[i]>1024.0)
        {
         pRatio[i]=1024.0;
        }
        cost_function[ite] = cost_function[ite] + (pFW[i] - pHist[i])*(pFW[i] - pHist[i]);
        //if(pFW[i] != pFW[i]) printf("--This is NaN at LOR=%d, with pHist=%f.\n", i, pHist[i]);
    }
    //printf("hiscount = %d, forcount = %d\n",hiscount,forcount);
}

void backward(float* pRatio, float* pRecon, float* pNorm, int t)
{
 int rx_plane = x_plane + 60;
 int index;
 float value;
 //x_plane=-1;
 TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));

 int i,j;

 FILE *inTAB,*inINCRE,*inSparse,*intSparse;

 for(i=0;i<4;i++)
 {
  for(j=0;j<4;j++)
  {
   char name[500];
   
   sprintf(kernel_name[i*4+j],"y%dz%d",j,i);

   sprintf(name, "../data_tof/sparse_table_plane_%d_%s_TOF",x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  
   int* tSparse= (int*) calloc(pTable[nccol_k].bias,sizeof(int)); 
   sprintf(name, "../data_tof/sparse_incre_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../data_tof/sparse_matrix_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   inSparse = fopen(name,"rb");
   fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);

   sprintf(name, "../data_tof/sparse_tof_plane_%d_%s_TOF", x_plane, kernel_name[i*4+j]);
   intSparse = fopen(name,"rb");
   fread(tSparse,sizeof(int),pTable[nccol_k].bias,intSparse);  

   int iz,jy,kcz,kcy;
   int LOR_initial=pTable[nccol_k].bias;
   for(iz=0;iz<ncz;iz++)
   {
   
    int shift_z = (iz*4-ncz*2)/4;
    for(jy=0;jy<ncy;jy++)
    {
   
      int shift_y= (jy*4-ncy*2)/4; 

     for(kcz=0;kcz<ncz;kcz++)
     {
      int cz1=kcz+52-shift_z;
      for(kcy=0;kcy<ncy;kcy++)
      {
       int cy1=kcy+68-shift_y;
       if(nboundary(cy1,cz1))
       {
        int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
        int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
        int c2total=pTable[cz1*ncy_k+cy1+1].bias;
        int ktotal;
        for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
        {
         int cy2=pIncre[ktotal]._y+cy2min;
         int cz2=pIncre[ktotal]._z+cz2min;
         if(boundary(cy2-68,cz2-52))
         {
	  if(tSparse[ktotal]==t)
	  {
          	pRecon[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pRatio[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]*pSparse[ktotal]/(2880*720);
          	pNorm[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pSparse[ktotal]/(2880*720);
	  }
         }
        }
       }
      }
     }
    }
   }
   fclose(inTAB);
   fclose(inINCRE);
   fclose(inSparse);
   fclose(intSparse);
   free(pIncre);
   free(pSparse);
   free(tSparse);
  }
 }
 free(pTable);
}

void update(float* pBwd,float* pRecon,float* pNorm,float* tvNorm,int ite)
{
 int ix,i;
 float beta;
 beta = 0;
 for(ix=0;ix<120;ix++)
 {
  for(i=0;i<nimgy*nimgz;i++)
  {
   if(pNorm[ix*nimgz*nimgy+i]+beta*tvNorm[ix*nimgz*nimgy+i])
   {
      pBwd[ix*nimgz*nimgy+i]=pBwd[ix*nimgz*nimgy+i]*pRecon[ix*nimgz*nimgy+i]/(pNorm[ix*nimgz*nimgy+i]+beta*tvNorm[ix*nimgz*nimgy+i]);
  
   }
  }
 }
}

int main()
{
 config();
 float* pHist= (float*) calloc(nccol*ncrow,sizeof(float));
 float* pFwd= (float*) calloc(nccol*ncrow,sizeof(float));
 float* pRatio= (float*) calloc(nccol*ncrow,sizeof(float));
 float* phantom= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pBwd= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pRecon= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pNorm= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* tvNorm= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* cost_function= (float*) calloc(100,sizeof(float));

 int i,ix,iy,iz;
 
 for(ix=0;ix<120;ix++)
 {
 	for(i=0;i<nimgy*nimgz;i++)
 	{
  		phantom[ix*nimgz*nimgy+i]=1.0;
		pBwd[ix*nimgz*nimgy+i]=1.0;
	}
 }


 memset(pFwd,0, ncrow*nccol * sizeof(float));
 memset(cost_function,0, 100*sizeof(float));

 
 FILE *iHist;
 FILE *ipBwd;
 //ipBwd=fopen("GPU_DHAPET_backward_TOF_ite5_.DAT","rb");
 ////ipBwd=fopen("GPU_DHAPET_backward_TOF_denoise_ite6_.DAT","rb");
 //fread(pBwd,sizeof(float),nimgx*nimgy*nimgz,ipBwd);
 //fclose(ipBwd);

 
/////////////////////////////////start em////////////////////////////////////////////
 int ite,tof,t;
 float T[2];
 char out_name[300];
 char norm_name[300];
 char tv_name[300];
 for(ite=0;ite<100;ite++)
 {
	T[0]=s_time();
        memset(pRecon,0,nimgx*nimgy*nimgz*sizeof(float));
  	memset(pNorm,0,nimgx*nimgy*nimgz*sizeof(float));
        
	for(tof=-1;tof<1;tof++)
	{
		t=abs(tof);
		char name[300];
		//sprintf(name, "../histogram/histogram_tof_%d.dat",tof);
                sprintf(name, "../histogram/histogram_attenuation_newnew.dat");
		iHist=fopen(name,"rb");
		memset(pHist,0,nccol*ncrow*sizeof(float));
 		fread(pHist,sizeof(float),nccol*ncrow,iHist);
 		fclose(iHist);
        	memset(pFwd,0,nccol*ncrow*sizeof(float));
		
  		for(x_plane=-1;x_plane<0;x_plane++)
  		{
  			forward(pFwd,pBwd,t,tof);
		}
  		ratio(pFwd,pHist,pRatio,ite, cost_function);
		for(x_plane=-1;x_plane<0;x_plane++)
		{
  			backward(pRatio,pRecon,pNorm,t);
		}
  		
	}
	printf("cost_function[%d] = %f\n",ite,sqrt(cost_function[ite]));
        update(pBwd,pRecon,pNorm,tvNorm,ite);
        memset(tvNorm,0,nimgx*nimgy*nimgz*sizeof(float));
        //tv_gradient(tvNorm,pBwd);
  	T[1]=s_time()-T[0];
  	printf("time: %g sec \n",T[1]);

        //if((ite+1)%10==0)
        //{
          sprintf(out_name,"attenuation_new_ite%d.DAT",ite+1);
	
	 
          FILE *fp2;
	 
	 
          if((fp2=fopen(out_name,"wb"))==NULL)
          {
	     printf("Can not open file\n");
 	     exit(0);
          }
          fwrite(pBwd,sizeof(float)*nimgx*nimgy*nimgz,1,fp2);
	  fclose(fp2);

	  /*sprintf(norm_name,"GPU_DHAPET_Norm_TV_ite%d_.DAT",ite+1);
	  FILE *ft;
          if((ft=fopen(norm_name,"wb"))==NULL)
          {
	     printf("Can not open file\n");
 	     exit(0);
          }
	  fwrite(pNorm,sizeof(float)*nimgx*nimgy*nimgz,1,ft);
	  fclose(ft);

	  sprintf(tv_name,"GPU_DHAPET_tvNorm_TV_ite%d_.DAT",ite+1);
	  FILE *fk;
          if((fk=fopen(tv_name,"wb"))==NULL)
          {
	     printf("Can not open file\n");
 	     exit(0);
          }
 	  fwrite(tvNorm,sizeof(float)*nimgx*nimgy*nimgz,1,fk);
          fclose(fk);*/
	 
	  
       //}
 }
/////////////////////////////////stop em///////////////////////////////////////////// 

 

 

 free(phantom); 
 free(pFwd);
 free(pBwd);
 free(pHist);
 free(pRatio);
 free(pRecon);
 free(pNorm);
 free(tvNorm);
}
