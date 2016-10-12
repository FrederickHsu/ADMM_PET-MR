#include <math.h>
#include <string.h>
#include "psf_analytic.h"

#define SWAP(a,b) memcpy(&temp, &a, sizeof(TRACING_CRY)); memcpy(&a, &b,sizeof(TRACING_CRY)); memcpy(&b, &temp, sizeof(TRACING_CRY)); 
#define MEMCPY(a,b) memcpy(&a, &b, sizeof(TRACING_CRY));
#define MINSERT 7
#define NSTACK 50

int istack[NSTACK];

void quicksort_tracing_cry(int n, TRACING_CRY arr[]) //arr[0...n-1]
{ //index of original code from Numerical Recipes is from 1.
    int i,ir = n-1,j,k,l = 0;
    int jstack = 0;

    TRACING_CRY a, temp;
 
    for(;;)
    {
        if(ir - l < MINSERT)
        {   
            for(j = l + 1; j <= ir; j++)
            {
                MEMCPY(a,arr[j])
                for(i = j - 1;i >= l; i--)
                {
                  if(arr[i].alpha <= a.alpha) break;
                  MEMCPY(arr[i+1],arr[i]);
                }
                MEMCPY(arr[i+1],a)
            }
            if(!jstack)
            {
                return;
            }
            ir = istack[jstack--];
            l  = istack[jstack--];
        }
        else
        {   
            k = (l + ir)>>1;
            SWAP(arr[k],arr[l+1])
            if(arr[l].alpha > arr[ir].alpha )
            { 
                SWAP(arr[l],arr[ir])
            } //SWAP is replacement of several lines, so should use "{...}" in if statement
            if(arr[l+1].alpha > arr[ir].alpha )
            {
                SWAP(arr[l+1],arr[ir])
            }
            if(arr[l].alpha > arr[l+1].alpha)
            {
                SWAP(arr[l],arr[l+1])
            }
            i = l + 1;
            j = ir;
            MEMCPY(a,arr[l+1])
            for(;;)
            {
                do { i++; } while(arr[i].alpha < a.alpha);
                do { j--; } while(arr[j].alpha > a.alpha);
                if(j < i) break;
                SWAP(arr[i], arr[j])
            }
            MEMCPY(arr[l+1],arr[j])
            MEMCPY(arr[j],a)
            jstack += 2;
            //Push pointers to larger subarray on stack, process smaller subarray immediately
            if(ir - i + 1 >= j)
            {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j - 1;
            }
            else
            {
                istack[jstack]   = j - 1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}

int merge_with_pitch(int n, TRACING_CRY arr[])
{
    int k, ns = 0;
    TRACING_CRY *pSTACK;

    for(k = 1; k < n; k++) ns = (arr[k].alpha == arr[k-1].alpha) ? (ns + 1) : ns;
    pSTACK = arr;
    for(k = 1; k < n; k++)
    {
        if(arr[k].alpha != pSTACK->alpha)  
        {
            memcpy(++pSTACK, arr + k, sizeof(TRACING_CRY));
        }
        else
        {
            if(arr[k].io != pSTACK->io)
            {
                pSTACK->io = 'o';
            }
        }
    }

    return n - ns;
}

int merge_without_pitch(int n, TRACING_CRY arr[])
{
    int k, ns = 0;
    TRACING_CRY *pSTACK;

    for(k = 1; k < n; k++) ns = (arr[k].alpha == arr[k-1].alpha) ? (ns + 1) : ns;
    pSTACK = arr;
    for(k = 1; k < n; k++)
    {
        if(arr[k].alpha != pSTACK->alpha)  
        {
            memcpy(++pSTACK, arr + k, sizeof(TRACING_CRY));
        }
    }

    return n - ns;
}
