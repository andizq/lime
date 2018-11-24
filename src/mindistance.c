/*
#include <math.h>
#include "constants.h"
#include "readdata.h"
#include "mindistance.h"
*/
#include "lime.h"

int standard_min_gp(double x,
		    double y,
		    double z, 
		    int Npoints, struct grid *gp){
  double mindist,distsq;
  mindist=1e5*PC;
  unsigned int i,j,k,ireal=0,ind=-1;

  for( i = 0; i < Npoints; i++){
    ireal = ID_picked[i];
    distsq = sqrt(  (x-gp[i].x[0])*(x-gp[i].x[0])	\
		  + (y-gp[i].x[1])*(y-gp[i].x[1])	\
		  + (z-gp[i].x[2])*(z-gp[i].x[2]));
		       
    if (distsq<mindist){
      mindist=distsq;
      ind=ireal;
    }
  }
  //printf("id picked %d out of %d, gp[ind].xyz : %.2lf, %.2lf, %.2lf, mindist %.4lf\n",ind,Npoints,gp[ind].x[0]/PC,gp[ind].x[1]/PC,gp[ind].x[2]/PC, mindist/PC);
  return ind;
}

int standard_min(double x, double *xm,
		 double y, double *ym,
		 double z, double *zm){
  double mindist,distsq;
  mindist=1e5*PC;
  unsigned int i,i0,j,k,ind=-1;

  i0 = 0;
  //printf("Starting point %d %d\n",i0);
  for( i = i0; i < Ndata; i++){
    distsq = sqrt(  (x-xm[i])*(x-xm[i])	\
		  + (y-ym[i])*(y-ym[i])	\
		  + (z-zm[i])*(z-zm[i]));
		       
    if (distsq<mindist){
      mindist=distsq;
      ind=i;
    }
  }
  //printf("mindist %.2lf, %.2lf, %.2lf\n",x,y,z);
  //printf("mindist %.2lf, id %d\n",mindist,ind);
  return ind;
}

int index_min(double u, double *um, int Nu){

  double mindist,distu;
  mindist=1e5*PC; //Huge initial value
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

  //extern unsigned short Nx,Ny,Nz; /* Already defined at readdata.c and declared at lime.h via readdata.h */
  int i,j,k,Num;

  i = index_min(x, xm, Nx);
  j = index_min(y, ym, Ny);
  k = index_min(z, zm, Nz);
  Num = i*Ny*Nz + j*Nz + k;
  
  return Num;
}
