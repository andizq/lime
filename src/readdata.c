#include "readdata.h"
#include "kdtree.h"

double **defcols2read(double **data){ //unsigned short *cols) {

  unsigned short noo, count = 0, count_abund = 0, i_abund = 0, foo = 4242;
  //  double **data;

  FILE *check_columns = fopen("header.dat", "r");
  fscanf(check_columns,"%hd",&foo);
  if (foo != 0){
    printf("ERROR (sf3dmodels input): Missing ids column, it MUST be provided by the user.\n");
    exit(1);
  }
  count++;
  for(noo = 1; noo < SF3D_max_cols; noo++){    
    fscanf(check_columns,"%hd",&foo);
    if (foo == SF3D_abundance) count_abund++;
    if (foo == 4242) break;
    count++;
  }
  //data = (double **)malloc(count * sizeof(double *));

  if (sf3d->cols != count){
    printf("ERROR (sf3dmodels input): Number of columns to be read ('header.dat') differs from the number of columns written on 'datatab.dat' ('npoints.dat')\n");
    printf("on 'npoints.dat': %d, on 'header.dat': %d\n",sf3d->cols,count);
    exit(1);
  }

  sf3d->abundance = (double **)malloc(count_abund * sizeof(double *));

  FILE *alloc_columns = fopen("header.dat", "r");

  for(noo = 0; noo < sf3d->cols; noo++ ){    

    fscanf(alloc_columns,"%hd",&foo);

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
    if (foo == SF3D_temp_gas){
      sf3d->temp_gas = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->temp_gas;
    }
    if (foo == SF3D_temp_dust){
      sf3d->temp_dust = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->temp_dust;
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
      sf3d->abundance[i_abund] = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->abundance[i_abund];
      i_abund++;
    }
    if (foo == SF3D_gtdratio){
      sf3d->gtdratio = malloc (sizeof(double) * Ndata);
      data[noo] = sf3d->gtdratio;
    }
    
  }
  
  fclose(check_columns);
  fclose(alloc_columns);
  
  /* Returning: */
  //*cols = count;
  
  //return data;

}

void readDatatab2() {

  unsigned int noo, i, *id;
  unsigned short j, sf3dcols;
  double **data;
  
  printf("*** sf3dmodels: %s, fixed_grid: %s\n",sf3dmodels ? "True" : "False", fixed_grid ? "True" : "False");
  printf("*** Looking for sf3dmodels input...\n");

  if(!fixed_grid){

    FILE *gridsize = fopen("npoints.dat", "r");
    fscanf(gridsize,"%hd %hd %hd %hd %d",&sf3dcols,&Nx,&Ny,&Nz,&Ndata);
    
    sf3d = malloc(sizeof(struct sf3d_data));  
    data = (double **)malloc(sf3dcols * sizeof(double *));
    //data = defcols2read();
    sf3d->cols = sf3dcols;
    defcols2read(data);
    id = sf3d->id;
    
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
    printf("   (Grid info. --> Nx,Ny,Nz,N: %hd %hd %hd %d)\n",Nx,Ny,Nz,Ndata);
    
    for (i = 0; i < Ndata; i++) {
      fscanf(fp, "%d", &id[i]);
      for (j = 1; j < sf3dcols; j++) {
	fscanf(fp, "%lf", &data[j][i]);
      }
    }
  
    printf("*** The data was succesfully read from 'x.dat' 'y.dat' 'z.dat' 'npoints.dat' 'datatab.dat'\n");

    fclose(fx);
    fclose(fy);
    fclose(fz);
    fclose(gridsize);
    fclose(fp);
    
  }else{
  
    FILE *gridsize = fopen("npoints.dat", "r");
    fscanf(gridsize,"%hd %hd %hd %hd %d",&sf3dcols,&Nx,&Ny,&Nz,&Ndata);
    
    sf3d = malloc(sizeof(struct sf3d_data));  
    //data = defcols2read(&cols);
    data = (double **)malloc(sf3dcols * sizeof(double *));
    sf3d->cols = sf3dcols;
    defcols2read(data);
    id = sf3d->id;

    FILE *fp  = fopen("datatab.dat", "r");
    
    printf("*** Found it. Reading the data...\n");

    for (i = 0; i < Ndata; i++) {
      fscanf(fp, "%d", &id[i]);
      for (j = 1; j < sf3dcols; j++) {
	fscanf(fp, "%lf", &data[j][i]);
      }
    }

    /*
    //new: next block
    //Insert the points into the KDTree object
    
    kd = kd_create(3); //new  
    unsigned int *datak, id_dat;
    datak = malloc (sizeof(unsigned int) * Ndata);
    printf("*** Inserting points into kdtree...\n");
    unsigned int nox, noy, noz;
    for( id_dat = 0; id_dat < Ndata; id_dat++ ){
      datak[id_dat] = id_dat;
      //printf("%d\n",datak[id_dat]);
      assert(kd_insert3(kd, 
			sf3d->x[id_dat], sf3d->y[id_dat], sf3d->z[id_dat], 
			&datak[id_dat]) == 0);
    }
    */


    printf("*** The data was succesfully read from 'npoints.dat' 'datatab.dat'\n");
    
    fclose(gridsize);
    fclose(fp);
    //printf("%f,%d,%d\n",sf3d->x[Ndata-1],cols,Ndata);
  }
  
}


void readDatatab() {

  unsigned int noo;
  
  printf("*** Looking for sf3dmodels input...\n");
  FILE *gridsize = fopen("npoints.dat", "r");
  fscanf(gridsize,"%hd %hd %hd %d",&Nx,&Ny,&Nz,&Ndata);
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

  //new: next block
  //Insert the points into the KDTree object
  /*
  kd = kd_create(3); //new  
  unsigned int *datak, id_dat=0;
  datak = malloc (sizeof(unsigned int) * Ndata);
  printf("*** Inserting points into kdtree...\n");
  unsigned int nox, noy, noz;
  for( nox = 0; nox < Nx; nox++ )
    for( noy = 0; noy < Ny; noy++ )
      for( noz = 0; noz < Nz; noz++ ){
	datak[id_dat] = id_dat;
	//printf("%d\n",datak[id_dat]);
	assert(kd_insert3(kd, xm[nox], ym[noy], zm[noz], &datak[id_dat]) == 0);
	id_dat++;
      }
  */

  //
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
  printf("   (Grid info. --> Nx,Ny,Nz,N: %hd %hd %hd %d)\n",Nx,Ny,Nz,Ndata);

  for( noo = 0; noo < Ndata; noo++ ){

    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf",
	   &ID[noo],&DENS[noo],&TEMP[noo],&VEL_x[noo],&VEL_y[noo],&VEL_z[noo],&ABUND[noo],&GTD[noo]);
    
  }
  
  printf("*** The data was succesfully read from 'x.dat' 'y.dat' 'z.dat' 'npoints.dat' 'datatab.dat'\n");

  fclose(fx);
  fclose(fy);
  fclose(fz);
  fclose(fp);
  
}

void freeDatatab() {
  free(DENS);
}
