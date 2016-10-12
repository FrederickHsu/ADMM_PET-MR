#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#define _GNU_SOURCE
#include <getopt.h>

#include "global_var_dec.h"
#include "psf_analytic.h"

struct option longopts [] = {
    {"crystal",1,NULL,'c'},
    {"debug",0,NULL,'d'},
    {"kernel",1,NULL,'k'},
    {"object",1,NULL,'o'},
    {"Spacing",1,NULL,'S'},
    {"xplane",1,NULL,'x'},
    {"DOI",0,NULL,'D'},
    {"inplane",1,NULL,'I'},
    {"pixel",1,NULL,'p'},
    {"grid",1,NULL,'g'},
    {"HRRT",1,NULL,'H'},
    {"help",0,NULL,'h'},
    {0,0,0,0}};

void help()
{
    printf("     -h,                        --help                               display this help.\n");
    printf("     -c <Num1:Num2:Num3:Num4>,  --crystal=Num1:Num2:Num3:Num4        give crystal pair, y1:z1:y2:z2.\n");
    printf("     -d,                        --debug                              debug mode.\n");
    printf("     -k <Num>,                  --kernel=Num                         kernel 0 ~ 2.\n");
    printf("     -S <Num>,                  --Spacing=Num                        Spacing size, mm.\n");
    printf("     -o <Str>,                  --object=Str                         description of saved files name.\n");
    printf("     -D,                        --DOI                                point spread function with depth of interation.\n");
    printf("     -x <Num>,                  --xplane=Num                         one specific xplane be computed for kernel.\n");
    printf("     -I <Num>,                  --inplane=Num                        numbers of inplanes be used in the image space.\n");
    printf("     -p <Num>,                  --pixel=Num                          how is the ratio crystal size to the size of pixel in one in-plane.\n");
    printf("     -H <NumCrslY:NumCrslZ>,    --HRRT=NumCrslY:NumCrslZ             the number crystals in each HRRT head.\n");
    printf("     -g <Num>,                  --grid=Num                           numbers of grids used in each voxel for computing PSF.\n");
}

void config()
{
    int i, j, k; // loop control variables
    int nkernel_triangle_stack;
    nkernel_triangle_stack = ((n_crsl2pxl % 2) == 0) ? (n_crsl2pxl / 2) : ((n_crsl2pxl + 1) / 2);
    /////////////////////////////////////////////////////////////////////////////////////////////
    //nkernel = (nkernel_triangle_stack + 1) * nkernel_triangle_stack / 2;
    nkernel=4;
    /////////////////////////////////////////////////////////////////////////////////////////////
    if(nkernel > 100)
    {
        printf("Too big number in -p <Num>, try small one. \n");
        exit(0);
    }

    nimgy = n_crsl2pxl * ncy;
    nimgz = n_crsl2pxl * ncz;    
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
    printf("  the object data will be saved with description: %s\n",description_name);
    py_edge_low = -ncy_k * cw / 2.;
    py_edge_up  =  ncy_k * cw / 2.;
    pz_edge_low = -ncz_k * cw / 2.;
    pz_edge_up  =  ncz_k * cw / 2.;

    if(x_plane >= nimgx / 2 || x_plane <= -nimgx / 2)
    {
        printf("Invalid number. In -x <Num>, the Num must be between %d and -%d.\n", nimgx/2-1, nimgx/2-1);
        help();
        exit(0);
    }

    vx = spacing / nimgx;
    vy = cw / n_crsl2pxl;
    vz = cw / n_crsl2pxl;

    kvx_edge_low = vx * (x_plane - 0.5);
    kvx_edge_up  = vx * (x_plane + 0.5);

    for(i = 0, j = 0; i < nkernel_triangle_stack; i++)
    {
        for(k = 0; k < nkernel_triangle_stack - i; k++, j++)
        {
            kvy_edge_low[j] = vy * (k + i);
            kvy_edge_up[j] = vy * (k + i + 1);
            kvz_edge_low[j] = vz * i;
            kvz_edge_up[j] = vz * (i + 1);
            ///////////////////////////////////////////////////////////////////////////////////////////////
            //sprintf(kernel_name[j],"y%dz%d", k+i, i);
            ///////////////////////////////////////////////////////////////////////////////////////////////
            //printf("%s\n", kernel_name[j]);

        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //for(j = 0; j < nkernel; j++)
    for(i = 0; i < nkernel; i++)
    {
        for(j = 0; j < nkernel; j++)
        {
         kernel_corner[i*nkernel+j][0]=j*vy;//+vy/2.0;
         kernel_corner[i*nkernel+j][1]=i*vz;//+vz/2.0;
         sprintf(kernel_name[i*nkernel+j],"y%dz%d",j,i);
        //kernel_corner[j][0] = kvy_edge_low[j];
        //kernel_corner[j][1] = kvz_edge_low[j];
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    pitch_half = pitch / 2;

}

void data_DIR()
{
    DIR *dp;
    if((dp = opendir("data")) == NULL )
    {
        printf("No directory: %s, and let me create it now.\n", "data");
        mkdir("data", S_IRUSR | S_IWUSR | S_IXUSR);
    }
    else
    {
        closedir(dp);
    }

}

int boundary(int cy, int cz)
{
    if( cy < 208 && cy >= 0 &&
        cz < 208 && cz >= 0    )
        return 1;
    else return 0;
}

int main(int argc, char *argv[])
{
    int opt;
    char tag_DOI = 1, tag_debug = 0;
    int icy1=ncy/2, icz1=ncz/2, icy2=ncy/2, icz2=ncz/2, ikernel=0;    

    while((opt = getopt_long(argc, argv, "c:dk:o:S:x:I:p:g:H:hD", longopts,NULL)) != -1)
    {
        switch(opt)
        {
            case 'h':
                help();
                exit(0);
                break;
            case 'o':
                strcpy(description_name,optarg);
                break;
            case 'k':
                ikernel = (float)atoi(optarg);
                if(ikernel < 0 && ikernel > 2)
                {
                    printf("Invalid number. In -k <Num>, Num is 0 ~ 2.\n");//, ncy - 1, ncz - 1);
                    exit(0);
                }
                break;
            case 'd':
                tag_debug = 1;
                break;
            case 'S':
                spacing = (float)atoi(optarg);
                break;
            case 'x':
                x_plane = atoi(optarg);
                break;
            case 'c':
                sscanf(optarg,"%d:%d:%d:%d",&icy1, &icz1, &icy2, &icz2);
                if(!(boundary(icy1, icz1) && boundary(icy2, icz2)))
                {
                    printf("Invalid number. In -d <Num1:Num2:Num3:Num4>, Num1 and Num3 are 0 ~ %d, Num2 and Num4 are 0 ~ %d.\n", ncy_k - 1, ncz_k - 1);
                    exit(0);
                }
                break;
            case 'D':
                tag_DOI = 1;
                break;
            case 'I':
                nimgx = atoi(optarg);
                break;
            case 'g':
                ngrid = atoi(optarg);
                break;
            case 'p':
                n_crsl2pxl = atoi(optarg);
                break;
            case 'H':
                sscanf(optarg,"%d:%d",&ncy, &ncz);
                break;
            case ':':
                printf("option needs a value\n");
                break;
            case '?':
                printf("unknown option: %c\n", optopt);
                break;
        }
    }

    data_DIR();
    config();

    /*if(tag_debug)
    {
        if(!tag_DOI)
        {
            debug_single_pair(icy1, icz1, icy2, icz2, ikernel);
            //debug_CONE(icy1, icz1, ikernel);
        }
        exit(0);
    }*/

    if(tag_DOI)
    {
        for(x_plane=-59;x_plane<60;x_plane++)
        {
        ray_tracing_panel();  // with DOI
        }
    }
    /*else
    {
        tracing_kvoxel();  // without DOI
    }*/

    exit(0);
}
