#include "readdata.h"

double **defcols2read(unsigned short *cols) {

  unsigned short noo, count = 0, foo = 4242;
  double **data;

  FILE *check_columns = fopen("npoints_test.dat", "r");
  for(noo = 0; noo < SF3D_max_cols; noo++ ){    
    fscanf(check_columns,"%d",&foo);
    if (noo == 0 && foo != 0){
      printf("ERROR (sf3dmodels input): Missing ids column, it MUST be provided by the user.\n");
      exit(1);
    }
    if (foo == 4242) break;
    count++;
  }
  data = (double **)malloc(count * sizeof(double *));

  FILE *alloc_columns = fopen("npoints_test.dat", "r");
  //fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);
  for(noo = 0; noo < count; noo++ ){    

    fscanf(alloc_columns,"%d",&foo);

    if (foo == SF3D_id){
      sf3d->id = malloc (sizeof(unsigned int) * Ndata);      
    }
    if (foo == SF3D_x){
      sf3d->x = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->x;
    }
    if (foo == SF3D_y){
      sf3d->y = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->y;
    }
    if (foo == SF3D_z){
      sf3d->z = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->z;
    }
    if (foo == SF3D_dens_H2){
      sf3d->dens_H2 = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_H2;
    }
    if (foo == SF3D_dens_p_H2){
      sf3d->dens_p_H2 = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_p_H2;
    }
    if (foo == SF3D_dens_o_H2){
      sf3d->dens_o_H2 = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_o_H2;
    }
    if (foo == SF3D_dens_e){
      sf3d->dens_e = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_e;
    }
    if (foo == SF3D_dens_H){
      sf3d->dens_H = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_H;
    }
    if (foo == SF3D_dens_He){
      sf3d->dens_He = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_He;
    }
    if (foo == SF3D_dens_Hplus){
      sf3d->dens_Hplus = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->dens_Hplus;
    }
    if (foo == SF3D_temperature){
      sf3d->temperature = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->temperature;
    }
    if (foo == SF3D_tdust){
      sf3d->tdust = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->tdust;
    }
    if (foo == SF3D_vel_x){
      sf3d->vel_x = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->vel_x;
    }
    if (foo == SF3D_vel_y){
      sf3d->vel_y = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->vel_y;
    }
    if (foo == SF3D_vel_z){
      sf3d->vel_z = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->vel_z;
    }
    if (foo == SF3D_abundance){
      sf3d->abundance = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->abundance;
    }
    if (foo == SF3D_gtdratio){
      sf3d->gtdratio = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->gtdratio;
    }

  }
  
  fclose(check_columns);
  fclose(alloc_columns);
  
  /* Returning: */
  *cols = count;
  return data;

}

void readDatatab2() {

  unsigned int noo, i, *id;
  unsigned short j, cols = 0;
  double **data;

  printf("sf3dmodels: %d, fixed_grid: %d\n",sf3dmodels,fixed_grid);
  printf("*** Looking for sf3dmodels input...\n");

  if(!fixed_grid){

  FILE *gridsize = fopen("npoints.dat", "r");
  fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);

  sf3d = malloc(sizeof(struct sf3d_data));  
  data = defcols2read(&cols);
  id = sf3d->id;
  sf3d->cols = cols;

  xm = malloc (sizeof(double) * Nx);
  FILE *fx  = fopen("x.dat", "r");
  for( noo = 0; noo < Nx; noo++ ) fscanf(fx,"%lf",&xm[noo]);

  ym = malloc (sizeof(double) * Ny);
  FILE *fy  = fopen("y.dat", "r");
  for( noo = 0; noo < Ny; noo++ ) fscanf(fy,"%lf",&ym[noo]);

  zm = malloc (sizeof(double) * Nz);
  FILE *fz  = fopen("z.dat", "r");
  for( noo = 0; noo < Nz; noo++ ) fscanf(fz,"%lf",&zm[noo]);
  
  FILE *fp  = fopen("datatab.dat", "r");  
  
  printf("*** Found it. Reading the data...\n");
  printf("   (Grid info. --> Nx,Ny,Nz,N: %d %d %d %d)\n",Nx,Ny,Nz,Ndata);
  
  for (i = 0; i < Ndata; i++) {
    fscanf(fp, "%d", &id[i]);
    for (j = 1; j < cols; j++) {
      fscanf(fp, "%lf", &data[j][i]);
    }
  }
  
  printf("*** The data was read succesfully from 'x.dat' 'y.dat' 'z.dat' 'npoints.dat' 'datatab.dat'\n");

  fclose(fx);
  fclose(fy);
  fclose(fz);
  fclose(gridsize);
  fclose(fp);

  }else{
  
    FILE *gridsize = fopen("npoints.dat", "r");
    fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);
    
    sf3d = malloc(sizeof(struct sf3d_data));  
    data = defcols2read(&cols);
    id = sf3d->id;
    sf3d->cols = cols;
    
    xm = malloc (sizeof(double) * Nx);
    FILE *fx  = fopen("x.dat", "r");
    for( noo = 0; noo < Nx; noo++ ) fscanf(fx,"%lf",&xm[noo]);
    
    ym = malloc (sizeof(double) * Ny);
    FILE *fy  = fopen("y.dat", "r");
    for( noo = 0; noo < Ny; noo++ ) fscanf(fy,"%lf",&ym[noo]);
    
    zm = malloc (sizeof(double) * Nz);
    FILE *fz  = fopen("z.dat", "r");
    for( noo = 0; noo < Nz; noo++ ) fscanf(fz,"%lf",&zm[noo]);
    
    FILE *fp  = fopen("datatab.dat", "r");
    
    printf("*** Found it. Reading the data...\n");

    for (i = 0; i < Ndata; i++) {
      fscanf(fp, "%d", &id[i]);
      for (j = 1; j < cols; j++) {
	fscanf(fp, "%lf", &data[j][i]);
      }
    }

    printf("*** The data was read succesfully from 'npoints.dat' 'datatab.dat'\n");

    fclose(fx);
    fclose(fy);
    fclose(fz);
    fclose(gridsize);
    fclose(fp);
    //printf("%f,%d,%d\n",sf3d->x[Ndata-1],cols,Ndata);
  }

  
}


void readDatatab() {

  unsigned int noo;
  
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
