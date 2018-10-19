#include "readdata.h"

void defcols2read() {
  unsigned short noo;
  unsigned short foo = 4242;
  //const int bsz=1; char buf[bsz];
  struct DATA sf3d; 
    
  FILE *gridsize = fopen("npoints_test.dat", "r");
  //fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);
  for(noo = 0; noo < SF3D_max_cols; noo++ ){    
    fscanf(gridsize,"%d",&foo);
    if (noo == 0 && foo != 0){
      printf("ERROR (sf3dmodels input): Missing ids column, it MUST be provided by the user.\n");
      exit(1);
    }
    if (foo == SF3D_id){ 
      sf3d.id = malloc (sizeof(unsigned int) * Ndata);
      printf("allocating\n");
      printf("%s\n",&sf3d.id);
    }
    if (foo == 4242){
      printf("breaking\n");
      break;
    }
  }
}

void readDatatab() {

  unsigned int noo;

  //defcols2read();

  printf("*** Looking for sf3dmodels input...\n");
  FILE *gridsize = fopen("npoints.dat", "r");
  fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);
  xm = malloc (sizeof(double) * Nx);
  FILE *fx  = fopen("x.dat", "r");
  for( noo = 0; noo < Nx; noo++ ){

    fscanf(fx,"%lf",&xm[noo]);
   
  }
 
  /////////////////////////////////////////

  ym = malloc (sizeof(double) * Ny);
  FILE *fy  = fopen("y.dat", "r");
  for( noo = 0; noo < Ny; noo++ ){

    fscanf(fy,"%lf",&ym[noo]);
   
  }

  /////////////////////////////////////////

  zm = malloc (sizeof(double) * Nz);
  FILE *fz  = fopen("z.dat", "r");
  for( noo = 0; noo < Nz; noo++ ){

    fscanf(fz,"%lf",&zm[noo]);
    
  }

  /////////////////////////////////////////

  ID = malloc (sizeof(unsigned int) * Ndata);
  DENS = malloc (sizeof(double) * Ndata);
  TEMP = malloc (sizeof(double) * Ndata);
  VEL_x = malloc (sizeof(double) * Ndata);
  VEL_y = malloc (sizeof(double) * Ndata);
  VEL_z = malloc (sizeof(double) * Ndata);
  ABUND = malloc (sizeof(double) * Ndata);
  GTD = malloc (sizeof(double) * Ndata);

  FILE *fp  = fopen("datatab.dat", "r");

  printf("*** Found it. Reading the data...\n");
  printf("   (Grid info. --> Nx,Ny,Nz,N: %d %d %d %d)\n",Nx,Ny,Nz,Ndata);

  for( noo = 0; noo < Ndata; noo++ ){

    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf",
	   &ID[noo],&DENS[noo],&TEMP[noo],&VEL_x[noo],&VEL_y[noo],&VEL_z[noo],&ABUND[noo],&GTD[noo]);
    
  }
  
  printf("*** The data was read succesfully from 'x.dat' 'y.dat' 'z.dat' 'npoints.dat' 'datatab.dat'\n");

  fclose(fx);
  fclose(fy);
  fclose(fz);
  fclose(fp);
  
}

void freeDatatab() {
  free(DENS);
}
