#ifndef READDATA_H
#define READDATA_H

#include <stdio.h>
#include <stdlib.h>
#include <collparts.h>
#include <assert.h>

#define SF3D_id             0
#define SF3D_x              1
#define SF3D_y              2
#define SF3D_z              3
#define SF3D_dens_H2        CP_H2 + 3     //4
#define SF3D_dens_p_H2      CP_p_H2 + 3   //5
#define SF3D_dens_o_H2      CP_o_H2 + 3   //6
#define SF3D_dens_e         CP_e + 3      //7
#define SF3D_dens_H         CP_H + 3      //8
#define SF3D_dens_He        CP_He + 3     //9
#define SF3D_dens_Hplus     CP_Hplus + 3  //10
#define SF3D_temperature    11
#define SF3D_tdust          12
#define SF3D_vel_x          13
#define SF3D_vel_y          14
#define SF3D_vel_z          15
#define SF3D_abundance      16
#define SF3D_gtdratio       17
#define SF3D_max_cols       18
//add temp_gas, temp_dust instead of just temperature

struct sf3d_data{
  unsigned int *id;
  double *x, *y, *z;
  double *dens_H2, *dens_p_H2, *dens_o_H2, *dens_e, *dens_H, *dens_He, *dens_Hplus; 
  double *temperature, *tdust;
  double *vel_x;
  double *vel_y;
  double *vel_z;
  double *abundance;
  double *gtdratio;
  unsigned short cols;
};


_Bool sf3dmodels, fixed_grid;
void *kd; //new
unsigned short Nx, Ny, Nz;
unsigned int Ndata, *ID, *ID_picked, *ids_fixed;
//extern unsigned int *ID_picked;
double *xm, *ym, *zm;
double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;
struct sf3d_data *sf3d;
void readDatatab(), readDatatab2(), freeDatatab();

#endif
