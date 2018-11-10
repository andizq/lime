#include <math.h>
#include "constants.h"
#include "readdata.h"
#include "mindistance.h"

double index_min(double u, double *um, int Nu){

  double mindist,distu;

  mindist=1000000*AU;

  //printf("Hey from index_min() at mindistance.c\n",mindist);
  /*if (mindist < par->radius){
    printf("%f\n",mindist);
  }
  */

  unsigned int i,j;
  for( i = 0; i < Nu; i++){
    distu = fabs(u-um[i]);

    if (distu<mindist){

      mindist=distu;
      j=i;
    }
  }
  return j;
}

int find_id_min(double x, double *xm, 
		double y, double *ym, 
		double z, double *zm){

  //extern unsigned short Nx,Ny,Nz; /* Already defined at readdata.c and thus at lime.h via readdata.h */
  int i,j,k,Num;

  i = index_min(x, xm, Nx);
  j = index_min(y, ym, Ny);
  k = index_min(z, zm, Nz);
  Num = i*Ny*Nz + j*Nz + k;
  
  return Num;
}
