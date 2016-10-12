#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#define _GNU_SOURCE
#include <getopt.h>
#define OUTFILE1 "DHAPET_original_forward.DAT"
#define OUTFILE2 "DHAPET_original_backward.DAT"
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

void forward(float* pFwd,float* phantom)
{
 int jlayer1=0;
 int jlayer2=0;

 int fx_plane = x_plane +60;
 TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));

 int i,j;

 FILE *inTAB,*inINCRE,*inSparse;


 for(i=0;i<4;i++)
 {
  for(j=0;j<4;j++)
  {
   char name[500];
   //char kernel_name[100][10];
   sprintf(kernel_name[i*4+j],"y%dz%d",j,i);

   sprintf(name, "../system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
   inSparse = fopen(name,"rb");
   fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse); 

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
  //if(c2total>0) printf("c2total=%d\n",c2total);
 // if(pTable[cz1*ncy_k/4+cy1].bias>0) printf("ktotal=%d\n",pTable[cz1*ncy_k/4+cy1].bias);
        int ktotal;
        for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
        {//printf("123\n");
         int cy2=pIncre[ktotal]._y+cy2min;
         int cz2=pIncre[ktotal]._z+cz2min;
//printf("cy2= %d cz2= %d ",cy2,cz2);

         if(boundary(cy2-68,cz2-52))
         {
          pFwd[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal]*phantom[fx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j];///(float)LOR_initial;
         // printf(" p: %f \n",phantom[59*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]);
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
   free(pIncre);
   free(pSparse);
  }
 }
 free(pTable);

} 

void ratio(float *pFW, float *pHist, float *pRatio)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {
        if(pHist[i] >= 0)
        {
            if(pFW[i] > 0)
            {
                pRatio[i] =  pFW[i]-pHist[i];
		if (pRatio[i]<0)
		{
			pRatio[i]=-pRatio[i];
		} 
		pRatio[i]=2*pRatio[i];
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
        //if(pFW[i] != pFW[i]) printf("--This is NaN at LOR=%d, with pHist=%f.\n", i, pHist[i]);
    }
}

void backward(float* pRatio, float* pRecon, float* pNorm)
{
 int rx_plane = x_plane +60;
 int jlayer1=0;
 int jlayer2=0;
 TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));

 int i,j;

 FILE *inTAB,*inINCRE,*inSparse;

 for(i=0;i<4;i++)
 {
  for(j=0;j<4;j++)
  {
   char name[500];
   //char kernel_name[100][10];
   sprintf(kernel_name[i*4+j],"y%dz%d",j,i);

   sprintf(name, "../system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
   inSparse = fopen(name,"rb");
   fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse); 

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
          pRecon[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pRatio[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]*pSparse[ktotal];
          pNorm[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pSparse[ktotal];
     
 
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
   free(pIncre);
   free(pSparse);
  }
 }
 free(pTable);
}



void proxtv(float*pRecon,float* pBwd,float* sol,float* x_n)

{
	float * r=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * s=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * dx=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * dy=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * p=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * q=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * pold=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * qold=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * temp_r=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * temp_s=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * sum_rs=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * weights=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * r2=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * s2=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * Y_n=(float*)calloc(nimgy*nimgz,sizeof(float));
	float tn=1,tn1;
	int i,j,k;
	float weight_x=1,weight_y=1,gamma=0.5,iteration=5,mt=1,t,told=1;

/////////////////u_n-gamma*f2.grad(u_n)///////////////////////////////
	for (i=0;i<nimgy*nimgz;i++)
	{	
		Y_n[i]=pBwd[i]-gamma*pRecon[i];
/////////////////initial////////////////////////////////////////////
		r[i]=0;
		s[i]=0;
		pold[i]=r[i];
		qold[i]=s[i];
	}


/////////////////////proxTV iteration//////////////////////////////
	for (j=0;j<iteration;j++)
	{
		for (i=0;i<nimgy*nimgz;i++)
		{
			if (i%nimgz==0)
			s2[i]=s[i];
			else if (i%nimgz==nimgz-1)
			s2[i]=-s[i-1];
			else
			s2[i]=s[i]-s[i-1];
			
			s2[i]=s2[i]*weight_x;

		}	
	
		for (i=0;i<nimgy;i++)
		{
		for (k=0;k<nimgz;k++)
		{
		if (i==0)
		r2[i*nimgy+k]=r[i*nimgy+k];
		else if (i==nimgy-1)
		r2[i*nimgy+k]=r[i*(nimgy-1)+k];	
		else 
		r2[i*nimgy+k]=r[i*nimgy+k]-r[i*(nimgy-1)+k];
		r2[i]=r2[i]*weight_y;
	
		}
		}
/////////////     current solution//////////////
		for (i=0;i<nimgy*nimgz;i++)
		{
		sol[i]=Y_n[i]-gamma*(r2[i]+s2[i]);
	
		}
////////////////////////////////////////////////////////////////////////////////////////////
		//update divergence vectors and project
//////////////////start gradient_op/////////////////////////////////		
		for (i=0;i<nimgy*nimgz;i++)
		{

		if (i%nimgz==nimgz-1)
		dy[i]=0;
		else
		dy[i]=sol[i+1]-sol[i];

		}

		for (i=0;i<nimgy;i++)
		{
		for (k=0;k<nimgz;k++)
		{
		if (i==nimgy-1)
		dx[i*nimgy+k]=0;
		else
		dx[i*nimgy+k]=sol[(i+1)*nimgy+k]-sol[i*nimgy+k];
		
		}
		}
//////////////////////////////////////////////////////////////////////////////////////////// 
		for (i=0;i<nimgy*nimgz;i++)

		{
		r[i]=r[i]-(1/(8*gamma)/(mt*mt)*dx[i]);
		s[i]=s[i]-(1/(8*gamma)/(mt*mt)*dy[i]);
		temp_r[i]=abs(r[i]);
		temp_s[i]=abs(s[i]);
		sum_rs[i]=temp_r[i]*temp_r[i]+temp_s[i]*temp_s[i];
		sum_rs[i]=sqrt(sum_rs[i]);
		if (sum_rs[i]<1)
			weights[i]=1;
		else 
		weights[i]=sum_rs[i];
		p[i]=r[i]/weights[i];
		q[i]=s[i]/weights[i];


		}
		
		//FISTA update
		t=(1+sqrt(4*told*told))/2;

		for (i=0;i<nimgy*nimgz;i++)
		{	
		r[i]=p[i]+((told-1)/t)*(p[i]-pold[i]);
		pold[i]=p[i];
		s[i]=q[i]+((told-1)/t)*(q[i]-qold[i]);
		qold[i]=q[i];


		}

		told=t;

	
	}



//	////////FISTA algorithm/////////
	tn1=(1+sqrt(1+4*tn*tn))/2;
	for (i=0;i<nimgy*nimgz;i++)
	{
	pBwd[i]=x_n[i]+((tn-1)/tn1)*(x_n[i]-sol[i]);

//	///////update//////////
	sol[i]=x_n[i];

	}	
	tn=tn1;

	free(r);
	free(s);	
	free(dx);
	free(dy);
	free(p);
	free(q);
	free(qold);
	free(pold);
	free(temp_r);
	free(temp_s);
	free(sum_rs);
	free(weights);
	free(r2);
	free(s2);
	free(Y_n);


}	














void update(float* pBwd,float* pRecon,float* pNorm)
{
 int ix,i;
 for(ix=0;ix<120;ix++)
 {
  for(i=0;i<nimgy*nimgz;i++)
  {
   if(pNorm[ix*nimgz*nimgy+i])
   {
    pBwd[ix*nimgz*nimgy+i]=pRecon[ix*nimgz*nimgy+i]/pNorm[ix*nimgz*nimgy+i];
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
 float* sol= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* x_n= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pRecon= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pNorm= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));

float* phantom_yz= (float*) calloc(nimgy*nimgz,sizeof(float));


 int i,ix,iy,iz;
 for(ix=0;ix<120;ix++)
 {
 	for(i=0;i<nimgy*nimgz;i++)
 	{
  		pBwd[ix*nimgz*nimgy+i]=0.0;
		x_n[ix*nimgz*nimgy+i]=0.0;
		sol[ix*nimgz*nimgy+i]=0.0;
	}
 }


 memset(pFwd,0, ncrow*nccol * sizeof(float));
 printf("stop 01\n");
 
 FILE *iHist;
 iHist=fopen("histogram.dat","rb");
 fread(pHist,sizeof(float),nccol*ncrow,iHist);
 fclose(iHist);
  printf("stop 02\n");
/////////////////////////////////start em////////////////////////////////////////////
 int ite;
 float T[2];
 char out_name[300];
 for(ite=0;ite<2;ite++)
 {
  T[0]=s_time();
  memset(pFwd,0,nccol*ncrow*sizeof(float));
  memset(pRecon,0,nimgx*nimgy*nimgz*sizeof(float));
  memset(pNorm,0,nimgx*nimgy*nimgz*sizeof(float));
  for(x_plane=-5;x_plane<0;x_plane++)
  {
	printf("start forward: %d\n",x_plane);  	
	forward(pFwd,pBwd);
  }
	
 printf("start ratio\n");
  ratio(pFwd,pHist,pRatio);

  for(x_plane=-5;x_plane<0;x_plane++)
  {
	//printf("start backward: %d\n",x_plane);  	  	
	backward(pRatio,pRecon,pNorm);
	proxtv(pRecon,pBwd,sol,x_n);
  }
  //update(pBwd,pRecon,pNorm);
  T[1]=s_time()-T[0];
  printf("time: %g sec \n",T[1]);

  if((ite+1)%1==0)
  {
    sprintf(out_name,"GPU_DHAPET_backward_ite%d_.DAT",ite+1);
    FILE *fp2;
    if((fp2=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
 	exit(0);
    }
    fwrite(pBwd,sizeof(float)*nimgx*nimgy*nimgz,1,fp2);
    fclose(fp2);
   }
 }
/////////////////////////////////stop em///////////////////////////////////////////// 

 

 //free(phantom_yz); 
 free(x_n);
 free(sol);
 free(phantom); 
 free(pFwd);
 free(pBwd);
 free(pHist);
 free(pRatio);
 free(pRecon);
 free(pNorm);
}
