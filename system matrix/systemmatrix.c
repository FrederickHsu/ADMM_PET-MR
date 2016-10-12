#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <time.h>
#include "global_var_dec.h"
#include "psf_analytic.h"

#define PI 3.1415926



//////////////////////////////////////

float s_time()
{
 return (float)clock()/CLOCKS_PER_SEC;
}


int is_hit_panel(float y, float z)
{
	if(y < py_edge_up - pitch_half && y > py_edge_low && z < pz_edge_up - pitch_half && z > pz_edge_low) return 1;
	else return 0;
}



	
void coincidence_response(float *yz_p1,float *yz_p2, float x, float y, float z, float *pCOIN, int *counter)
{
	int CRY_ID_upy, CRY_ID_upz, CRY_ID_lowy, CRY_ID_lowz, CIU, CIL ;
	CRY_ID_upy = CRY_ID_upz = CRY_ID_lowy = CRY_ID_lowz = CIU = CIL =0;
	float y_p1,y_p2,z_p1,z_p2; 
	
	y_p1=yz_p1[0];
	z_p1=yz_p1[1];
	y_p2=yz_p2[0];
	z_p2=yz_p2[1];

	CRY_ID_lowy = floor((y_p1-py_edge_low)/cw);
	CRY_ID_lowz = floor((z_p1-pz_edge_low)/cw);
	CRY_ID_upy = floor((y_p2-py_edge_low)/cw);
	CRY_ID_upz = floor((z_p2-pz_edge_low)/cw);
	
	
	counter[0]=counter[0]+1;
	pCOIN[(CRY_ID_lowz * ncy_k + CRY_ID_lowy) * nccol_k + CRY_ID_upz * ncy_k + CRY_ID_upy]+=1.0;
	
		

	return;
	
}





int sparse_matrix_wr(TABLE_ZIP* pTable, INCREMENT_ZIP* pIncre, float* pSparse, float* TOP_cone_resp, FILE *outf_incre, FILE *outf_matrix)
{
    
    int total;
    int cymin, cymax, czmin, czmax;
    int irow, jcol;
    int cy2, cz2, y, z;
    int time;
    cymin = ncy_k - 1;
    cymax = 0;
    czmin = ncz_k - 1;
    czmax = 0;

    for(jcol = 0; jcol < nccol_k; jcol++)
    {
	//time = TOP_cone_resp[jcol];
	
        if(TOP_cone_resp[jcol] > 0)
        {
	    //printf("TOP_cone_resp[jcol]=%d\n",time);
            cz2 = jcol / ncy_k; // ncy_k = ncz_k
            cy2 = jcol % ncz_k;
            czmin = (czmin > cz2) ? cz2 : czmin;
            czmax = (czmax < cz2) ? cz2 : czmax;
            cymin = (cymin > cy2) ? cy2 : cymin;
            cymax = (cymax < cy2) ? cy2 : cymax;
        }
    }

    pTable->cymin = cymin;
    pTable->czmin = czmin;

    for(z = czmin, total = 0; z <= czmax; z++)
    {
        for(y = cymin, jcol = z * ncy_k + cymin; y <= cymax; y++, jcol++)
        {
            if(TOP_cone_resp[jcol] > 0)
            {
                pIncre[total]._y = y - cymin;
                pIncre[total]._z = z - czmin;
                pSparse[total] = TOP_cone_resp[jcol];
                //dSparse[total] = DOI_cone_resp[jcol];
                total++; 
            }
        }
    }
    
    int nw_incre, nw_matrix;
    

    if(total)
    {
	printf("total=%d\n",total);
        nw_incre=fwrite(pIncre, sizeof(INCREMENT_ZIP), total, outf_incre);
	//printf("total=%d, nw=%d\n",total, nw_incre);
        nw_matrix=fwrite(pSparse, sizeof(float), total, outf_matrix);
        //printf("total=%d, incre=%d, matrix=%d\n",total,nw_incre,nw_matrix);
    }

    return total;

}




	
void ray_tracing_panel()
{
	float y_in_p, z_in_p;
	float yz_p1[2], yz_p2[2]; // recode the position that trace comes in each layer of each panel
    	//yz_p*[0][0] = y_in_p ;
    	//yz_p*[0][1] = z_in_p ; * = 1 or 2
	
    	// The above: Start point and end point  between ray and panel.
	float PHI, THETA;
    	float cPHI, sPHI, cTHETA, sTHETA;
    	int iPHI, iTHETA, nTHETA, nPHI;
    	nPHI=720;
    	float fly1_L[2], fly2_L[2];
	float L_bw_x;
	float *pCOIN = (float*) calloc(ncrow_k * nccol_k, sizeof(float));

    	
//--------------------------------------------------------------------
	FILE *outf, *outf_incre, *outf_matrix, *outf_table;
    	int irow, bias;
    	char name[500];

	float* TOP_cone_resp = (float*) calloc(nccol_k, sizeof(float));	
 	TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(ncrow_k +1, sizeof(TABLE_ZIP));
    	INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(nccol_k, sizeof(INCREMENT_ZIP));
    	float* pSparse = (float*) calloc(nccol_k, sizeof(float));  
//--------------------------------------------------------------------
	

	float y0, z0, x0;
    	float y, z, x;
    	float lgridy = vy / ngrid, lgridz = vz / ngrid, lgridx = vx / ngrid;
    	int igx, igy, igz, k;
    	int ikernel;
    	int count;
	int counter[4];

	//x0 = (x_plane - 0.5) * vx + lgridx / 2;
    	x0= x_plane * vx + lgridx / 2;
	float T[2];
        T[0]=s_time();

	
	for(ikernel = 0; ikernel < nkernel*nkernel; ikernel++)
	{

	    counter[0]=counter[1]=counter[2]=counter[3]=0;
	    printf("\n  computing the analytic PSF of kernel %s\n", kernel_name[ikernel]);

            sprintf(name,"data/sparse_incre_plane_%d_%s_", x_plane, kernel_name
[ikernel]); outf_incre = fopen(name, "wb");
            sprintf(name,"data/sparse_matrix_plane_%d_%s_", x_plane, kernel_name[ikernel]); outf_matrix = fopen(name, "wb");
            sprintf(name,"data/sparse_table_plane_%d_%s_", x_plane, kernel_name[ikernel]); outf_table = fopen(name, "wb");
	   y0 = kernel_corner[ikernel][0] + lgridy / 2;
           z0 = kernel_corner[ikernel][1] + lgridz / 2;
	
	   memset(pCOIN, 0, sizeof(float) * ncrow_k * nccol_k);
	   bias = 0;

	   ///////////////
	   float T[2];
           T[0]=s_time();
	   //////////////
	
	   for(igx = 0, x = x0; igx < ngrid; igx++, x += lgridx)
           { 
              for(igz = 0, z = z0; igz < ngrid; igz++, z += lgridz)
              {
                 for(igy = 0, y = y0; igy < ngrid; igy++, y += lgridy)
                 {
                    //printf("pass%d, \n",igz*ngrid+igy);
                    for(iPHI = 0, count = 0; iPHI < nPHI; iPHI++)
                    {
		       nTHETA = (iPHI>0)?2880:0;

                       PHI = iPHI * PI / (2 * nPHI); // PHI is sampled at 0.125 degree.
                       //float radius = iPHI * ( ncz * vz ) / nPHI;
                       //float height = (x >= 0)?abs(-spacing/2.0 - x):(spacing/2.0 - x);
                       //float d_radius_1=(x>=0)?(ncz*vz)*((spacing-height)/height)/nPHI:ncz*vz/nPHI;
                       //float d_radius_2=(x>=0)?ncz*vz/nPHI:(ncz*vz)*((spacing-height)/height)/nPHI;
                       //PHI = atan(radius/height);

                       cPHI = cos(PHI);
                       sPHI = sin(PHI);
                       fly1_L[0] = fabs((spacing / 2 - x) / cPHI);
                       //fly1_L[1] = fabs((spacing / 2 + 2 * cl - x) / cPHI);
                       fly2_L[0] = fabs((-spacing / 2 - x) / cPHI);
                       //fly2_L[1] = fabs((-spacing / 2 - 2 * cl - x) / cPHI);
		       //float circle= (spacing / 2.0)*sin(PHI)*2.0*PI;
                       //float circle=(radius>0)?radius*2.0*PI:0.1*PI;
                       //nTHETA = (circle>0.1*PI)?(circle/(0.005*PI)):1;
                       //nTHETA = circle/(0.01*PI);
                       //nTHETA=floor(nTHETA);
                        	 
                       L_bw_x = 2 * cl / cPHI;
		       for(iTHETA = 0; iTHETA < nTHETA; iTHETA++)
                       {
                           THETA = iTHETA * PI *2 / nTHETA; // THETA is sampled at 0.125 degree.
                           cTHETA = cos(THETA);
                           sTHETA = sin(THETA);

                           //float d_circle_1=(x>=0)?circle*((spacing-height)/height)/nTHETA:circle/nTHETA;
                           //float d_circle_2=(x>=0)?circle/nTHETA:circle*((spacing-height)/height)/nTHETA;

                           y_in_p   = y + sPHI * cTHETA * fly1_L[0];
                           z_in_p   = z + sPHI * sTHETA * fly1_L[0];
                           yz_p1[0] = y_in_p ;
                           yz_p1[1] = z_in_p ;
                            		
			   y_in_p   = y - sPHI * cTHETA * fly2_L[0];
                           z_in_p   = z - sPHI * sTHETA * fly2_L[0];
                           yz_p2[0] = y_in_p ;
                           yz_p2[1] = z_in_p ;
                           if(is_hit_panel(yz_p1[0], yz_p1[1]) && is_hit_panel(yz_p2[0], yz_p2[1]))
			   {
				
			 	coincidence_response(yz_p1, yz_p2, x, y, z, pCOIN, counter);
				
			   }//end if
			}//end for iTHETA
		      }//end for iPHi
		   }//end for igy
		}//end for igz
	     }//end for igx

	    // printf("coincidence=%d,tof_0=%d,tof_1=%d,tof_2=%d\n",counter[0],counter[1],counter[2],counter[3]);
	     
	     for(irow = 0; irow < ncrow_k; irow++)
             {
            	 pTable[irow].bias = bias;
        	 memcpy(TOP_cone_resp, pCOIN + irow * nccol_k, nccol_k * sizeof(float));
            	 bias += sparse_matrix_wr(pTable + irow, pIncre, pSparse, TOP_cone_resp, outf_incre, outf_matrix);
             }
	     

	     pTable[irow].bias = bias;
             fwrite(pTable, sizeof(TABLE_ZIP), ncrow_k + 1, outf_table);

	     T[1]=s_time()-T[0];
             printf("time: %g sec \n",T[1]);
	     fclose(outf_incre);
             fclose(outf_matrix);
             fclose(outf_table);

	}//end for ikernel


	free(pCOIN);
    	free(TOP_cone_resp);
    	free(pTable);
    	free(pIncre);
    	free(pSparse);



}
