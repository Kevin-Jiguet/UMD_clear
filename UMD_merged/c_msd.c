#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void compute_msd(const double *pos, double *msd, const int hh, const int vv, const int ballistic, const int nitmax)
{

        for (int i=ballistic+hh ; i<nitmax ; i=i+hh){
		
            const double Xorigin = pos[3*i];
            const double Yorigin = pos[3*i+1];
            const double Zorigin = pos[3*i+2];

	    

            for (int j=0 ; j<nitmax ; j=j+vv){

                const double Xfuture = pos[3*(i+j)];
                const double Yfuture = pos[3*(i+j)+1];
                const double Zfuture = pos[3*(i+j)+2];

                msd[(j/vv)] = msd[(j/vv)] + pow((Xfuture-Xorigin),2) + pow((Yfuture-Yorigin),2) + pow((Zfuture-Zorigin),2);

            }
        }
}
