//
// Created by hsa on 20/12/16.
//

#ifndef SJI_SIMULATION_H
#define SJI_SIMULATION_H

#include "sjinc.h"

#define C50 ( 1.239407e+0f)
#define C51 (-1.105315e-1f)
#define C52 ( 2.496329e-2f)
#define C53 (-5.804879e-3f)
#define C54 ( 9.358680e-4f)

#define sjmsgfd2dn1(p, ix, iz) (C50*(p[ix][iz+1]-p[ix][iz-0]) + \
                                C51*(p[ix][iz+2]-p[ix][iz-1]) + \
                                C52*(p[ix][iz+3]-p[ix][iz-2]) + \
                                C53*(p[ix][iz+4]-p[ix][iz-3]) + \
                                C54*(p[ix][iz+5]-p[ix][iz-4]) )

#define sjmsgfd2dn2(p, ix, iz) (C50*(p[ix+1][iz]-p[ix-0][iz]) + \
                                C51*(p[ix+2][iz]-p[ix-1][iz]) + \
                                C52*(p[ix+3][iz]-p[ix-2][iz]) + \
                                C53*(p[ix+4][iz]-p[ix-3][iz]) + \
                                C54*(p[ix+5][iz]-p[ix-4][iz]) )


//! Survey
int sjssurvey_init(sjssurvey *ptr);

int sjssurvey_display(sjssurvey *ptr);

int sjssurvey_readis(sjssurvey *ptr, int is);

int sjssurvey_write(sjssurvey *ptr, int ifappend);

int sjssurvey_getparas(sjssurvey *ptr, int argc, char **argv);

//! Source
int sjssource_init(sjssource *ptr);

int sjssource_display(sjssource *ptr);

int sjssource_getparas(sjssource *ptr, int argc, char **argv);

//! model
int sjsgeo_init(sjsgeo *ptr);

int sjsgeo_display(sjsgeo *ptr);

int sjsgeo_getparas2d(sjsgeo *ptr, int argc, char **argv, char *info);

//! wave
int sjswave_init(sjswave *ptr);

int sjswave_display(sjswave *ptr);

int sjswave_getparas(sjswave *ptr, int argc, char **argv, char *info);

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

//! Two dimension acoustic simulation based on constant velocity-stress equation
void sjawsgfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

//! Two dimension acoustic reverse time simulation based on constant velocity-stress equation
void sjawrtsgfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

#endif //SJI_SIMULATION_H
