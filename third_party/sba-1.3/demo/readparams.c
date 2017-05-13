/* Loading of camera, 3D point & image projection parameters from disk files */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXSTRLEN  4096 /* 4K */

static char buf[MAXSTRLEN];
/* get rid of the rest of a line upto \n or EOF */
#define SKIP_LINE(f){                                                       \
  while(!feof(f))                                                           \
    if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

#include "eucsbademo.h"

/* returns the number of cameras defined in a camera parameters file.
 * Each line of the file corresponds to the parameters of a single camera
 */
static int readNcameras(FILE *fp)
{
int lineno, ncams, ch;

  lineno=ncams=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    SKIP_LINE(fp);
    ++lineno;
    if(ferror(fp)){
      fprintf(stderr, "readNcameras(): error reading input file, line %d\n", lineno);
      exit(1);
    }
    ++ncams;
  }

  return ncams;
}

/* reads into "params" the camera parameters defined in a camera parameters file.
 * "params" is assumed to point to sufficiently large memory.
 * Each line contains the parameters of a single camera, 7 parameters per camera are
 * assumed
 */
static void readCameraParams(FILE *fp, double *params)
{
const int cnp=7;
int lineno, n, ch;

  lineno=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);
    ++lineno;
    n=fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf\n", params, params+1, params+2, params+3, params+4, params+5, params+6); 
    if(n!=cnp){
      fprintf(stderr, "readCameraParams(): line %d contains %d parameters, expected %d!\n", lineno, n, cnp);
      exit(1);
    }
    if(ferror(fp)){
      fprintf(stderr, "readNcameras(): error reading input file, line %d\n", lineno);
      exit(1);
    }

    params+=cnp;
  }
}


/* determines the number of 3D points contained in a points parameter file as well as the
 * total number of their 2D image projections across all images. The file format is
 * X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ...
 */
static void readNpointsAndNprojections(FILE *fp, int *n3Dpts, int *nprojs)
{
int lineno, npts, nframes, ch;

  *n3Dpts=*nprojs=lineno=npts=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      ++lineno;
      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);
    ++lineno;
    fscanf(fp, "%*g%*g%*g%d", &nframes);
    if(ferror(fp)){
      fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d\n", lineno);
      exit(1);
    }
    SKIP_LINE(fp);
    *nprojs+=nframes;
    ++npts;
  }

  *n3Dpts=npts;
}


/* reads a points parameter file.
 * "params", "projs" & "vmask" are assumed preallocated, pointing to
 * memory blocks large enough to hold the parameters of 3D points, 
 * their projections in all images and the point visibility mask, respectively.
 * File format is X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ...
 */
static void readPointParamsAndProjections(FILE *fp, double *params, double *projs, char *vmask, int ncams)
{
int nframes, ch, lineno, ptno, frameno, n;
register int i;

  lineno=ptno=0;
  while(!feof(fp)){
    if((ch=fgetc(fp))=='#'){ /* skip comments */
      SKIP_LINE(fp);
      lineno++;

      continue;
    }

    if(feof(fp)) break;

    ungetc(ch, fp);

    fscanf(fp, "%lf%lf%lf", params, params+1, params+2); /* read in X Y Z */
    params+=3;

    fscanf(fp, "%d", &nframes); /* read in number of image projections */

    for(i=0; i<nframes; ++i){
      n=fscanf(fp, "%d%lf%lf", &frameno, projs, projs+1); /* read in image projection */
      if(n!=3){
        fprintf(stderr, "readPointParamsAndProjections(): error reading image projections from line %d [n=%d].\n"
                        "Line contains fewer than %d projections?\n", lineno+1, n, nframes);
        exit(1);
      }

      if(frameno>=ncams){
        fprintf(stderr, "readPointParamsAndProjections(): line %d contains an image projection for frame %d "
                        "but only %d cameras have been specified!\n", lineno+1, frameno, ncams);
        exit(1);
      }

      projs+=2;
      vmask[ptno*ncams+frameno]=1;
    }

    fscanf(fp, "\n"); // consume trailing newline

    lineno++;
    ptno++;
  }
}


/* combines the above routines to read the initial estimates of the motion + structure parameters from text files.
 * Also, it loads the projections of 3D points across images. The routine dynamically allocates the required amount
 * of memory (last 3 arguments).
 */
void readInitialSBAEstimate(char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs, double **motstruct, double **imgpts, char **vmask)
{
const int cnp=7, /* 4 rot params + 3 trans params */
          pnp=3, /* euclidean 3D points */
          mnp=2; /* img ponts are 2D */

FILE *fpc, *fpp;

  if((fpc=fopen(camsfname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", camsfname);
    exit(1);
  }

  if((fpp=fopen(ptsfname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", ptsfname);
    exit(1);
  }

  *ncams=readNcameras(fpc);
  readNpointsAndNprojections(fpp, n3Dpts, n2Dprojs);

  *motstruct=(double *)malloc((*ncams*cnp + *n3Dpts*pnp)*sizeof(double));
  if(*motstruct==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  *imgpts=(double *)malloc(*n2Dprojs*mnp*sizeof(double));
  if(*imgpts==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  *vmask=(char *)malloc(*n3Dpts * *ncams * sizeof(char));
  if(*vmask==NULL){
    fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
    exit(1);
  }
  memset(*vmask, 0, *n3Dpts * *ncams * sizeof(char)); /* clear vmask */


  /* prepare for re-reading files */
  rewind(fpc);
  rewind(fpp);

  readCameraParams(fpc, *motstruct);
  readPointParamsAndProjections(fpp, *motstruct+*ncams*cnp, *imgpts, *vmask, *ncams);

  fclose(fpc);
  fclose(fpp);
}

/* reads the 3x3 intrinsic calibration matrix contained in a file */
void readCalibParams(char *fname, double ical[9])
{
  FILE *fp;
  int i;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", fname);
    exit(1);
  }

  for(i=0; i<3; i++){
    fscanf(fp, "%lf%lf%lf\n", ical, ical+1, ical+2);
    ical+=3;
  }

  fclose(fp);
}

/* prints the initial estimates of the motion + structure parameters. It also prints the projections
 * of 3D points across images. For debugging purposes only.
 */
void printSBAData(double *motstruct, int ncams, int n3Dpts, double *imgpts, int n2Dprojs, char *vmask)
{
const int cnp=7, /* 4 rot params + 3 trans params */
          pnp=3, /* euclidean 3D points */
          mnp=2; /* img ponts are 2D */

register int i;

  printf("Motion parameters:\n");
  for(i=0; i<ncams*cnp; ++i)
    printf("%lf ", motstruct[i]);

  motstruct+=i;
  printf("\n\nStructure parameters:\n");
  for(i=0; i<n3Dpts*pnp; ++i)
    printf("%lf ", motstruct[i]);

   printf("\n\nImage projections:\n");
  for(i=0; i<n2Dprojs*mnp; ++i)
    printf("%lf ", imgpts[i]);

  printf("\n\nVisibility mask\n");
  for(i=0; i<ncams*n3Dpts; ++i)
    printf("%d%s", (int)vmask[i], ((i+1)%ncams)? " " : "\n");
  printf("\n");
}
