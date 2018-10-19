#ifndef READDATA_H
#define READDATA_H

#include <stdio.h>
#include <stdlib.h>
#include <collparts.h>

_Bool sf3dmodels;
unsigned short Nx, Ny, Nz;
unsigned int Ndata, *ID, *ID_picked;
double *xm, *ym, *zm;
double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;
void readDatatab(), freeDatatab();


#define SF3D_id             0
#define SF3D_dens_H2        CP_H2     //1
#define SF3D_dens_pH2       CP_p_H2   //2
#define SF3D_dens_oH2       CP_o_H2   //3
#define SF3D_dens_e         CP_e      //4
#define SF3D_dens_H         CP_H      //5
#define SF3D_dens_He        CP_He     //6
#define SF3D_dens_Hplus     CP_Hplus  //7
#define SF3D_temperature    8
#define SF3D_vel_x          9
#define SF3D_vel_y          10
#define SF3D_abundance      11
#define SF3D_gtdratio       12
#define SF3D_max_cols       13

struct DATA{
  unsigned int *id;
  double *dens_H2, *dens_pH2, *dens_oH2, *dens_e, *dens_H, *dens_He, *dens_Hplus; 
  double *temperature;
  double *vel_x;
  double *vel_y;
  double *vel_z;
  double *abundance;
  double *gtdratio; 
};

#endif
