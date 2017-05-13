#ifndef _SBA_DEMO_H_
#define _SBA_DEMO_H_

#define SBA_TERMINATION_THRESH  1E-12 //1E-15

#define SBA_MAX_REPROJ_ERROR    4.0 // max motion only reprojection error

#define BA_MOTSTRUCT            0
#define BA_MOT                  1
#define BA_STRUCT               2
#define BA_MOT_MOTSTRUCT        3

/* in imgproj.c */
extern void calcImgProj(double *a, double *qr, double *t, double *M, double *n);
extern void calcImgProjJacRTS(double *a, double *qr, double *t, double *M, double jacmRT[2][7], double jacmS[2][3]);

/* in readparams.c */
extern void readInitialSBAEstimate(char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs,
                                   double **motstruct, double **imgpts, char **vmask);
extern void readCalibParams(char *fname, double ical[9]);
extern void printSBAData(double *motstruct, int ncams, int n3Dpts, double *imgpts, int n2Dprojs, char *vmask);

/* in eucsbademo.c */
extern void sba_driver(char *camsfname, char *ptsfname, char *calibfname);

#endif /* _SBA_DEMO_H_ */
