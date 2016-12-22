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

#define B60 (-2.982778e+0f)
#define B61 ( 1.714286e+0f)
#define B62 (-2.678571e-1f)
#define B63 ( 5.291005e-2f)
#define B64 (-8.928571e-3f)
#define B65 ( 1.038961e-3f)
#define B66 (-6.012506e-5f)

#define B611  0.562500000000
#define B612 -0.112500000000
#define B613  0.012500000000
#define B622  0.022500000000
#define B623 -0.002500000000
#define B633  0.000277777778

#define sjmfd2dn1(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix][iz+1]+a[ix][iz-1]) + \
                            B62*(a[ix][iz+2]+a[ix][iz-2]) + \
                            B63*(a[ix][iz+3]+a[ix][iz-3]) + \
                            B64*(a[ix][iz+4]+a[ix][iz-4]) + \
                            B65*(a[ix][iz+5]+a[ix][iz-5]) + \
                            B66*(a[ix][iz+6]+a[ix][iz-6]) )

#define sjmfd2dn2(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix+1][iz]+a[ix-1][iz]) + \
                            B62*(a[ix+2][iz]+a[ix-2][iz]) + \
                            B63*(a[ix+3][iz]+a[ix-3][iz]) + \
                            B64*(a[ix+4][iz]+a[ix-4][iz]) + \
                            B65*(a[ix+5][iz]+a[ix-5][iz]) + \
                            B66*(a[ix+6][iz]+a[ix-6][iz]) )

#define sjmfd2dnc(a, ix, iz)( B611*(a[ix+1][iz+1]-a[ix-1][iz+1]-a[ix+1][iz-1]+a[ix-1][iz-1]) + \
                            B612*(a[ix+1][iz+2]-a[ix-1][iz+2]-a[ix+1][iz-2]+a[ix-1][iz-2]+a[ix+2][iz+1]-a[ix-2][iz+1]-a[ix+2][iz-1]+a[ix-2][iz-1]) + \
                            B613*(a[ix+1][iz+3]-a[ix-1][iz+3]-a[ix+1][iz-3]+a[ix-1][iz-3]+a[ix+3][iz+1]-a[ix-3][iz+1]-a[ix+3][iz-1]+a[ix-3][iz-1]) + \
                            B622*(a[ix+2][iz+2]-a[ix-2][iz+2]-a[ix+2][iz-2]+a[ix-2][iz-2]) + \
                            B623*(a[ix+2][iz+3]-a[ix-2][iz+3]-a[ix+2][iz-3]+a[ix-2][iz-3]+a[ix+3][iz+2]-a[ix-3][iz+2]-a[ix+3][iz-2]+a[ix-3][iz-2]) + \
                            B633*(a[ix+3][iz+3]-a[ix-3][iz+3]-a[ix+3][iz-3]+a[ix-3][iz-3]))

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

//! Two dimension acoustic simulation based on constant density equation
void sjawfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

//! Two dimension acoustic reverse time simulation based on constant density equation
void sjawrtfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

#endif //SJI_SIMULATION_H
