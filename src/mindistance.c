/*
#include <math.h>
#include "constants.h"
#include "readdata.h"
#include "mindistance.h"
*/
#include "lime.h"
#include "kdtree.h"

static double dist_sq( double *a1, double *a2, int dims ); //new

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
  
  /*
  int i,j,k,Num;
  
  i = index_min(x, xm, Nx);
  j = index_min(y, ym, Ny);
  k = index_min(z, zm, Nz);
  Num = i*Ny*Nz + j*Nz + k; //commented

  return Num;

  */

  int Num = 0;
  double pt[3] = { x, y, z };
  struct kdres *presults;
  char *data;
  unsigned int *pch, res_size = 0;
  double pos[3], dist;
  
  presults = kd_nearest( kd, pt );
  pch = (unsigned int*)kd_res_item( presults, pos );
  return *pch;
  
  /*
  presults = kd_nearest_range( kd, pt, radius_kd );
  while (1){
    res_size = kd_res_size(presults);
    if (res_size == 0) radius_kd = radius_kd * 1.5;
    else if (res_size > 10) radius_kd = radius_kd / 2.3;
    else break;
    presults = kd_nearest_range( kd, pt, radius_kd );
  }
  pch = (unsigned int*)kd_res_item( presults, pos );
  return *pch;
  */

  /*
  printf( "found %d results:\n", kd_res_size(presults) );
  dist = sqrt( dist_sq( pt, pos, 3 ) );
  printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data %d\n", 
	  pos[0], pos[1], pos[2], dist, *pch);  
  */
  
  /*
  presults = kd_nearest_range( kd, pt, radius );
  printf( "found %d results:\n", kd_res_size(presults) );

  while( !kd_res_end( presults ) ) {
    // get the data and position of the current result item 
    pch = (char*)kd_res_item( presults, pos );
    printf( "found %d results:\n", kd_res_size(presults) );
    // compute the distance of the current result from the pt 
    dist = sqrt( dist_sq( pt, pos, 3 ) );
printf( "found %.3f results:\n", dist );
// print out the retrieved data 
    printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data\n", 
    pos[0], pos[1], pos[2], dist);

    // go to the next entry 
    kd_res_next( presults );
  }
  */

  
}

static double dist_sq( double *a1, double *a2, int dims ) {
  double dist_sq = 0, diff;
  while( --dims >= 0 ) {
    diff = (a1[dims] - a2[dims]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}
