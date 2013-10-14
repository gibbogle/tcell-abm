#ifndef LIBPARA32_H
#define LIBPARA32_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void execute(int *,char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *,int *,int *);
void get_scene(int *, int *, int *, int *, int *, int *);
void get_summary(int *);
void get_profile_cd69(double *, double *, int *);
void get_profile_s1pr1(double *, double *, int *);
void get_profile_stim(double *, double *, int *);
void get_profile_stimrate(double *, double *, int *);
void get_profile_firstdccontacttime(double *, double *, int *);
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBPARA32_H
