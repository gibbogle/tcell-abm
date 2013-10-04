#ifndef TRANSFER_H
#define TRANSFER_H

#include <QMutex.h>

extern int showingVTK;

#define MAX_TC 1000000
#define MAX_DC 1000
#define MAX_BOND 50000

#define PROFILE_CD69 0
#define PROFILE_S1PR1 1
#define PROFILE_STIM 2
#define PROFILE_STIMRATE 3

extern int VTKbuffer[100];
extern int TC_list[5*MAX_TC];
extern int nTC_list;
extern int DC_list[5*MAX_DC];
extern int nDC_list;
extern int bond_list[2*MAX_BOND];
extern int nbond_list;
extern QMutex mutex1, mutex2;

extern int summaryData[100];
extern int NX, NY, NZ;
extern int nt_vtk;
extern bool leftb;

extern double *profile_x[10];
extern double *profile_y[10];
extern int profile_n[10];

#endif // TRANSFER_H
