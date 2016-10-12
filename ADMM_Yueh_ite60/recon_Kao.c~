#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
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

void norm_LOR(float* pNorm_LOR)
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

   sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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
          pNorm_LOR[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal];
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

void norm_voxel(float* pNorm_voxel)
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

   sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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
          
          pNorm_voxel[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pSparse[ktotal];
 
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


void forward(float* pFwd,float* phantom,float* pNorm_LOR)
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

   sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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
          pFwd[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal]*phantom[fx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]/pNorm_LOR[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)];///(float)LOR_initial;
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
         //   if(pFW[i] > 0)
            //{
                pRatio[i] =  abs(pFW[i]-pHist[i]);
		/*if (pRatio[i]<0)
		{
			pRatio[i]=-pRatio[i];
		} */
		pRatio[i]=2*pRatio[i];
           // }
           /* else
            {
                if(pHist[i] > 0)
                {
                    pRatio[i] = 0; // non_zero / zero
                }
                else
                { 
                    pRatio[i] = 1.0; // zero / zero
                }
            }*/
        }
       /*if(pRatio[i]>1024.0)
        {
         pRatio[i]=1024.0;
        }*/

    }

}

void norm(float *sol, float* norm2,int ite)
{
  int i;

  for(i=0;i<nimgx*nimgy*nimgz;i++)
  { 
     norm2[ite] = norm2[ite] + sol[i]*sol[i];
  }

}


void cost_function(float * pFW, float *pHist, float* norm2, float *costfunction, int ite, float gamma)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {

                costfunction[ite] =  costfunction[ite] + pow((pFW[i]-pHist[i]),2);		      
         
        
    }
    printf("first term[%d] = %f\n", ite, costfunction[ite]);
    printf("second term = %f\n",gamma*sqrt(norm2[ite]));
    costfunction[ite] =  costfunction[ite] + gamma *sqrt(norm2[ite]);
    printf("cost function[%d] = %f\n", ite, costfunction[ite]);
}

void backward(float* pRatio, float* pRecon, float* pNorm_voxel)
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

   sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inTAB = fopen(name,"rb");
   fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
   INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
   float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));  

   sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
   inINCRE = fopen(name,"rb");
   fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);

   sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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
          pRecon[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pRatio[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]*pSparse[ktotal]/pNorm_voxel[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j];     
 
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



void proxtv(float* x_n,float* u_n, float* r, float *s, float* r2, float* s2,float* dx, float* dy,float* pold, float* qold, float* p, float* q, float* temp_s, float* temp_r, float* sum_rs,float gamma, float* Y_n)

{
	float weights;
	
	float weight_x=1,weight_y=1,mt=1,t,told=1;
	int rx_plane = x_plane +60;
 	int i,k;
/////////////////////proxTV iteration//////////////////////////////
	
	//printf("%d",iteration);		
	for (i=0;i<nimgy*nimgz;i++)
		{
			if (i%nimgy==0)
			{	s2[i]=s[i];}
			else if (i%nimgy==nimgy-1)
			{	s2[i]=-s[i-1];}
			else
			{s2[i]=s[i]-s[i-1];}
			
			s2[i]=s2[i]*weight_x;

		}	
	
		for (i=0;i<nimgz;i++)
		{
		for (k=0;k<nimgy;k++)
		{
		if (i==0)
			{r2[i*nimgy+k]=r[i*nimgy+k];}
		else if (i==nimgz-1)
			{r2[i*nimgy+k]=-1*r[(i-1)*(nimgy)+k];}	
		else 
			{r2[i*nimgy+k]=r[i*nimgy+k]-r[(i-1)*nimgy+k];}
		
r2[i]=r2[i]*weight_y;
	
		}
		}
/////////////     current solution//////////////
		for (i=0;i<nimgy*nimgz;i++)
		{
		x_n[rx_plane*nimgy*nimgz+i]=Y_n[i]-gamma*(r2[i]+s2[i]);
		//sol[rx_plane*nimgy*nimgz+i]=Y_n[i]-gamma*(r2[i]+s2[i]);
		}
////////////////////////////////////////////////////////////////////////////////////////////
		//update divergence vectors and project
//////////////////start gradient_op/////////////////////////////////		
		for (i=0;i<nimgy*nimgz;i++)
		{

		if (i%nimgy==nimgy-1)
		dx[i]=0;
		else
		dx[i]=x_n[rx_plane*nimgy*nimgz+i+1]-x_n[rx_plane*nimgy*nimgz+i];

		}

		for (i=0;i<nimgz;i++)
		{
		for (k=0;k<nimgy;k++)
		{
		if (i==nimgz-1)
		dy[i*nimgy+k]=0;
		else
		dy[i*nimgy+k]=x_n[rx_plane*nimgy*nimgz+(i+1)*nimgy+k]-x_n[rx_plane*nimgy*nimgz+i*nimgy+k];
		
		}
		}
//////////////////////////////////////////////////////////////////////////////////////////// 
		float max=0;		
		float weight;
		for (i=0;i<nimgy*nimgz;i++)

		{
		r[i]=r[i]-(1/(8*gamma)/(mt*mt)*dy[i]);
		s[i]=s[i]-(1/(8*gamma)/(mt*mt)*dx[i]);
		temp_r[i]=abs(r[i]);
		temp_s[i]=abs(s[i]);
		sum_rs[i]=temp_r[i]*temp_r[i]+temp_s[i]*temp_s[i];
		sum_rs[i]=sqrt(sum_rs[i]);
			
		if (sum_rs[i]>=max)
		{
		max=sum_rs[i];
		
		}
		weights=max;
		if (max<1)
			weights=1;
		

		}

		for (i=0;i<nimgy*nimgz;i++)
		{

		//weights[i]=sum_rs[i];
		p[i]=r[i]/weights;
		q[i]=s[i]/weights;


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


	




	
	


/*
//	////////FISTA algorithm/////////
	tn1=(1+sqrt(1+4*tn*tn))/2;
	for (i=0;i<nimgy*nimgz;i++)
	{
	pBwd[rx_plane*nimgy*nimgz+i]=x_n[rx_plane*nimgy*nimgz+i]+((tn-1)/tn1)*(x_n[rx_plane*nimgy*nimgz+i]-sol[rx_plane*nimgy*nimgz+i]);

//	///////update//////////
	

	}	
	tn=tn1;
*/
	


}	














/*void update(float* sol,float* pNorm)
{
 int ix,i;
 for(ix=0;ix<120;ix++)
 {
  for(i=0;i<nimgy*nimgz;i++)
  {
   if(pNorm[ix*nimgz*nimgy+i])
   {
    //pBwd[ix*nimgz*nimgy+i]=pRecon[ix*nimgz*nimgy+i]/pNorm[ix*nimgz*nimgy+i];
    sol[ix*nimgz*nimgy+i]=sol[ix*nimgz*nimgy+i];
///pNorm[ix*nimgz*nimgy+i];
//)+0.5;	
   }
  }
 }
}*/

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
 float* u_n= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pRecon= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 
 float* pNorm_voxel= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 float* pNorm_LOR= (float*) calloc(nccol*ncrow,sizeof(float));
 float* costfunction= (float*) calloc(100,sizeof(float));
 float* norm2=(float*) calloc(100,sizeof(float));;

	float * Y_n=(float*)calloc(nimgy*nimgz,sizeof(float));
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
	//float * weights=(float*)calloc(nimgy*nimgz,sizeof(float));
	
	
	float * r2=(float*)calloc(nimgy*nimgz,sizeof(float));
	float * s2=(float*)calloc(nimgy*nimgz,sizeof(float));
	
     
	float tn=1,tn1;
	int i,j,k;
	float gamma=0.1;
	int iteration=100;
	int rx_plane = x_plane +60;



//float* phantom_yz= (float*) calloc(nimgy*nimgz,sizeof(float));
 int ix,iy,iz,kk,ss;

 memset(pFwd,0, ncrow*nccol * sizeof(float));
 memset(costfunction,0, 100*sizeof(float));
 memset(norm2,0, 100*sizeof(float));
 //printf("stop 01\n");
 
 FILE *iHist;
 iHist=fopen("../histogram/histogram_activity_newnew.dat","rb");//attenuation:histogram_att_1.dat;activity:histogram_att_211.dat
 fread(pHist,sizeof(float),nccol*ncrow,iHist);
 fclose(iHist);
// printf("stop 02\n");


 for(ix=0;ix<120;ix++)
 {
 	for(i=0;i<nimgy*nimgz;i++)
 	{
  		pBwd[ix*nimgz*nimgy+i]=0.0;
		x_n[ix*nimgz*nimgy+i]=pHist[ix*nimgz*nimgy+i];
               // x_n[ix*nimgz*nimgy+i]=0.0;
		sol[ix*nimgz*nimgy+i]=0.0;
		u_n[ix*nimgz*nimgy+i]=0.0;
	}
 }



/////////////////////////////////start prox////////////////////////////////////////////
 int ite;
 float T[2];
 char out_name[300];
 
 
  
  T[0]=s_time();
  memset(pFwd,0,nccol*ncrow*sizeof(float));
  memset(pRecon,0,nimgx*nimgy*nimgz*sizeof(float));
  memset(pNorm_voxel,0,nimgx*nimgy*nimgz*sizeof(float));
  memset(pNorm_LOR,0,nccol*ncrow*sizeof(float));
  
  for(x_plane=0;x_plane<3;x_plane++)
  {     
	x_plane=x_plane*(-20);
        norm_LOR(pNorm_LOR);
        norm_voxel(pNorm_voxel);
        x_plane=x_plane/(-20);
  }

  for(x_plane=0;x_plane<3;x_plane++)
  {     
	x_plane=x_plane*(-20);
	//printf("start forward: %d\n",x_plane);  	
	forward(pFwd,pBwd,pNorm_LOR);
        x_plane=x_plane/(-20);
  }
	
 //printf("start ratio\n");
  ratio(pFwd,pHist,pRatio);

  for(x_plane=0;x_plane<3;x_plane++)
  {     
	x_plane=x_plane*(-20);
	//printf("start backward: %d\n",x_plane);  	  	
	
	backward(pRatio,pRecon,pNorm_voxel);
        x_plane=x_plane/(-20);
  }


  for (ite=0;ite<iteration;ite++)
  {

/////////////////u_n-gamma*f2.grad(u_n)///////////////////////////////
     for(x_plane=0;x_plane<3;x_plane++)
     {	    
            x_plane=x_plane*(-20);
	    for (i=0;i<nimgy*nimgz;i++)
	    {	
		rx_plane = x_plane +60;
		Y_n[i]=pBwd[rx_plane*nimgy*nimgz+i]-gamma*pRecon[rx_plane*nimgy*nimgz+i];
/////////////////initial////////////////////////////////////////////
		r[i]=0;
		s[i]=0;
		pold[i]=r[i];
		qold[i]=s[i];
	    }

	    //printf("iteration=%d\n",ite+1);
	    proxtv(x_n,u_n,r,s,r2,s2,dx,dy,pold,qold,p,q,temp_s,temp_r,sum_rs,gamma,Y_n);
	    for(i=0;i<nimgy*nimgz;i++)
	    {
		sol[rx_plane*nimgy*nimgz+i]=x_n[rx_plane*nimgy*nimgz+i];
	    }
                
	    if(((ite+1)%5==0) || (ite<5))
 	    {   
		if(x_plane==-40)   
		{
        	   forward(pFwd,sol,pNorm_LOR);
        	   norm(sol,norm2,ite);
		   //printf("norm2=%f\n",norm2[ite]);
        	   cost_function(pFwd,pHist,norm2,costfunction,ite,gamma);
                //}
	        //printf("write file\n");

 	        //if((ite+1)%1==0)
 	        //{
	        //    if((ite+5)%5==0)
	        //    {
   		   sprintf(out_name,"GPU_DHAPET_backward_ite%d_activity_r01.DAT",ite+1);
   	 	   FILE *fp2;
   	 	   if((fp2=fopen(out_name,"wb"))==NULL)
    		   {
       	 	     printf("Can not open file\n");
 		     exit(0);
   		   }
   	 	   fwrite(sol,sizeof(float)*nimgx*nimgy*nimgz,1,fp2);
   	           fclose(fp2);
	         }
   	     }
	     x_plane=x_plane/(-20);
	 }
     }
     //update(pBwd,pRecon,pNorm);
     //update(sol,pNorm);
     T[1]=s_time()-T[0];
     printf("time: %g sec \n",T[1]);


/////////////////////////////////stop em///////////////////////////////////////////// 

 
 free(norm2);
 //free(phantom_yz); 
 free(x_n);
 free(u_n);
 free(sol);
 free(phantom); 
 free(pFwd);
 free(pBwd);
 free(pHist);
 free(pRatio);
 free(pRecon);
 free(pNorm_LOR);
 free(pNorm_voxel);
 free(costfunction);


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
	//free(weights);
	free(r2);
	free(s2);
	free(Y_n);








}
