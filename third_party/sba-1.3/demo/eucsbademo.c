/* Euclidean bundle adjustment demo using the sba package */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sba.h>
#include "eucsbademo.h"

#define CLOCKS_PER_MSEC (CLOCKS_PER_SEC/1000.0)


/* pointers to additional data, used for computed image projections and their jacobians */
struct globs_{
	double *intrcalib; /* intrinsic callibration matrix in row-major storage */

	double *ptparams; /* needed only when bundle adjusting for camera parameters only */
	double *camparams; /* needed only when bundle adjusting for structure parameters only */
} globs;


/* Routines to estimate the estimated measurement vector (i.e. "func") and
 * its sparse jacobian (i.e. "fjac") needed in BA. Code below makes use of the
 * routines calcImgProj() and calcImgProjJacRTS() which compute the predicted
 * projection & jacobian of a SINGLE 3D point (see imgproj.c).
 * In the terminology of TR-340, these routines compute Q and its jacobians A=dQ/da, B=dQ/db.
 * Notice also that what follows is two pairs of "func" and corresponding "fjac" routines.
 * The first is to be used in full (i.e. motion + structure) BA, the second in 
 * motion only BA.
 */


/**********************************************************************/
/* MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE SIMPLE DRIVERS */
/**********************************************************************/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given the parameter vectors aj and bi of camera j and point i, computes in xij the
 * predicted projection of point i on image j
 */
static void img_projRTS(int j, int i, double *aj, double *bi, double *xij, void *adata)
{
  double *Kcalib;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;

  calcImgProj(Kcalib, aj, aj+4, bi, xij); // 4 is the quaternion's length
}

/* Given the parameter vectors aj and bi of camera j and point i, computes in Aij, Bij the
 * jacobian of the predicted projection of point i on image j
 */
static void img_projRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata)
{
  double *Kcalib;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;

  calcImgProjJacRTS(Kcalib, aj, aj+4, bi, (double (*)[7])Aij, (double (*)[3])Bij); // 4 is the quaternion's length
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given the parameter vector aj of camera j, computes in xij the
 * predicted projection of point i on image j
 */
static void img_projRT(int j, int i, double *aj, double *xij, void *adata)
{
  const int pnp=3; /* euclidean 3D points */

  double *Kcalib, *ptparams;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  ptparams=gl->ptparams;

  calcImgProj(Kcalib, aj, aj+4, ptparams+i*pnp, xij); // 4 is the quaternion's length
}

/* Given the parameter vector aj of camera j, computes in Aij
 * the jacobian of the predicted projection of point i on image j
 */
static void img_projRT_jac(int j, int i, double *aj, double *Aij, void *adata)
{
  const int pnp=3; /* euclidean 3D points */

  double *Kcalib, *ptparams;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  ptparams=gl->ptparams;

  calcImgProjJacRTS(Kcalib, aj, aj+4, ptparams+i*pnp, (double (*)[7])Aij, NULL); // 4 is the quaternion's length
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given the parameter vector bi of point i, computes in xij the
 * predicted projection of point i on image j
 */
static void img_projS(int j, int i, double *bi, double *xij, void *adata)
{
  const int cnp=7;

  double *Kcalib, *camparams, *aj;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProj(Kcalib, aj, aj+4, bi, xij); // 4 is the quaternion's length
}

/* Given the parameter vector bi of point i, computes in Bij
 * the jacobian of the predicted projection of point i on image j
 */
static void img_projS_jac(int j, int i, double *bi, double *Bij, void *adata)
{
  const int cnp=7;

  double *Kcalib, *camparams, *aj;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProjJacRTS(Kcalib, aj, aj+4, bi, NULL, (double (*)[3])Bij); // 4 is the quaternion's length
}

/**********************************************************************/
/* MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS */
/**********************************************************************/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
static void img_projsRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *Kcalib;
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=pa+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProj(Kcalib, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
static void img_projsRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pa, *pb, *pqr, *pt, *ppt, *jaca, *jacb, *pA, *pB, *Kcalib;
  //int n;
  int m, nnz, Asz, Bsz, idx;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp;
  jaca=jac; jacb=jac+idxij->nnz*Asz;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=pa+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      idx=idxij->val[rcidxs[i]];
      pA=jaca + idx*Asz; // set pA to point to A_ij
      pB=jacb + idx*Bsz; // set pB to point to B_ij

      calcImgProjJacRTS(Kcalib, pqr, pt, ppt, (double (*)[7])pA, (double (*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB
    }
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given a parameter vector p made up of the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images.
 * The measurements are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T,
 * where hx_ij is the predicted projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
static void img_projsRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pqr, *pt, *ppt, *pmeas, *Kcalib, *ptparams;
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=p+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
	    ppt=ptparams + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProj(Kcalib, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the parameters of m cameras, compute in jac
 * the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm),
 * where A_ij=dx_ij/db_j (see HZ).
 * Notice that depending on idxij, some of the A_ij might be missing
 *
 */
static void img_projsRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pqr, *pt, *ppt, *pA, *Kcalib, *ptparams;
  //int n;
  int m, nnz, Asz, idx;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;
  Asz=mnp*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=p+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=ptparams + rcsubs[i]*pnp;
      idx=idxij->val[rcidxs[i]];
      pA=jac + idx*Asz; // set pA to point to A_ij

      calcImgProjJacRTS(Kcalib, pqr, pt, ppt, (double (*)[7])pA, (double (*)[3])NULL); // evaluate dQ/da in pA
    }
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
static void img_projsS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pqr, *pt, *ppt, *pmeas, *Kcalib, *camparams;
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=camparams+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProj(Kcalib, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the B_ij might be missing
 *
 */
static void img_projsS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  const int cnp=7, /* 4 rot params + 3 trans params */
            pnp=3, /* euclidean 3D points */
            mnp=2; /* img ponts are 2D */
  double *pqr, *pt, *ppt, *pB, *Kcalib, *camparams;
  //int n;
  int m, nnz, Bsz, idx;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kcalib=gl->intrcalib;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;
  Bsz=mnp*pnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=camparams+j*cnp;
    pt=pqr+4; // rot. quaternion has 4 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      idx=idxij->val[rcidxs[i]];
      pB=jac + idx*Bsz; // set pB to point to B_ij

      calcImgProjJacRTS(Kcalib, pqr, pt, ppt, (double (*)[7])NULL, (double (*)[3])pB); // evaluate dQ/da in pB
    }
  }
}


/* Driver for sba_xxx_levmar */
void sba_driver(char *camsfname, char *ptsfname, char *calibfname)
{
  double *motstruct, *motstruct_copy, *imgpts;
  double ical[9]; // intrinsic calibration params
  char *vmask;
  double opts[SBA_OPTSSZ], info[SBA_INFOSZ], phi;
  int howto, expert, analyticjac, n, verbose=0;
  int nframes, numpts3D, numprojs, nvars;
  const int cnp=7, pnp=3, mnp=2, nconstframes=1;

  static char *howtoname[]={"BA_MOTSTRUCT", "BA_MOT", "BA_STRUCT", "BA_MOT_MOTSTRUCT"};

  clock_t start_time, end_time;


  readInitialSBAEstimate(camsfname, ptsfname, &nframes, &numpts3D, &numprojs, &motstruct, &imgpts, &vmask);

  //printSBAData(motstruct, nframes, numpts3D, imgpts, numprojs, vmask);

  /* set up globs structure */
  readCalibParams(calibfname, ical); 
  globs.intrcalib=ical;
  globs.ptparams=NULL;
  globs.camparams=NULL;

  /* call sparse LM routine */
  opts[0]=SBA_INIT_MU; opts[1]=SBA_TERMINATION_THRESH; opts[2]=SBA_TERMINATION_THRESH;
  opts[3]=SBA_TERMINATION_THRESH;
  //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05

  /* Notice the various BA options demonstrated below */

  /* minimize motion & structure, motion only, o
   * motion and possibly motion & structure in a 2nd pass?
   */
  howto=BA_MOTSTRUCT;
  //howto=BA_MOT;
  //howto=BA_STRUCT;
  //howto=BA_MOT_MOTSTRUCT;

  /* simple or expert drivers? */
  //expert=0;
  expert=1;

  /* analytic or approximate jacobian? */
  //analyticjac=0;
  analyticjac=1;

  start_time=clock();
  switch(howto){
    case BA_MOTSTRUCT: /* BA for motion & structure */
      nvars=nframes*cnp+numpts3D*pnp;
      if(expert)
        n=sba_motstr_levmar_x(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, mnp,
                            img_projsRTS_x, (analyticjac)? img_projsRTS_jac_x : NULL, (void *)(&globs), 150, verbose, opts, info);
      else
        n=sba_motstr_levmar(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, mnp,
                            img_projRTS, (analyticjac)? img_projRTS_jac : NULL, (void *)(&globs), 150, verbose, opts, info);
    break;

    case BA_MOT: /* BA for motion only */
      globs.ptparams=motstruct+nframes*cnp;
      nvars=nframes*cnp;
      if(expert)
        n=sba_mot_levmar_x(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, imgpts, mnp,
                          img_projsRT_x, (analyticjac)? img_projsRT_jac_x : NULL, (void *)(&globs), 100, verbose, opts, info);
      else
        n=sba_mot_levmar(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, imgpts, mnp,
                          img_projRT, (analyticjac)? img_projRT_jac : NULL, (void *)(&globs), 100, verbose, opts, info);
    break;

    case BA_STRUCT: /* BA for structure only */
      globs.camparams=motstruct;
      nvars=numpts3D*pnp;
      if(expert)
        n=sba_str_levmar_x(numpts3D, nframes, vmask, motstruct+nframes*cnp, pnp, imgpts, mnp,
                          img_projsS_x, (analyticjac)? img_projsS_jac_x : NULL, (void *)(&globs), 100, verbose, opts, info);
      else
        n=sba_str_levmar(numpts3D, nframes, vmask, motstruct+nframes*cnp, pnp, imgpts, mnp,
                          img_projS, (analyticjac)? img_projS_jac : NULL, (void *)(&globs), 100, verbose, opts, info);
    break;

    case BA_MOT_MOTSTRUCT: /* BA for motion only; if error too large, then BA for motion & structure */
      if((motstruct_copy=(double *)malloc((nframes*cnp + numpts3D*pnp)*sizeof(double)))==NULL){
        fprintf(stderr, "memory allocation failed in readInitialSBAEstimate()\n");
        exit(1);
      }

      memcpy(motstruct_copy, motstruct, (nframes*cnp + numpts3D*pnp)*sizeof(double)); // save starting point for later use
      globs.ptparams=motstruct+nframes*cnp;
      nvars=nframes*cnp;

      if(expert)
        n=sba_mot_levmar_x(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, imgpts, mnp,
                          img_projsRT_x, (analyticjac)? img_projsRT_jac_x : NULL, (void *)(&globs), 100, verbose, opts, info);
      else
        n=sba_mot_levmar(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, imgpts, mnp,
                        img_projRT, (analyticjac)? img_projRT_jac : NULL, (void *)(&globs), 100, verbose, opts, info);

      if((phi=info[1]/numprojs)>SBA_MAX_REPROJ_ERROR){
        fflush(stdout); fprintf(stdout, "Refining structure (motion only error %g)...\n", phi); fflush(stdout);
        memcpy(motstruct, motstruct_copy, (nframes*cnp + numpts3D*pnp)*sizeof(double)); // reset starting point

        if(expert)
          n=sba_motstr_levmar_x(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, mnp,
                                img_projsRTS_x, (analyticjac)? img_projsRTS_jac_x : NULL, (void *)(&globs), 150, verbose, opts, info);
        else
          n=sba_motstr_levmar(numpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, mnp,
                              img_projRTS, (analyticjac)? img_projRTS_jac : NULL, (void *)(&globs), 150, verbose, opts, info);
      }
      free(motstruct_copy);

    break;

    default:
      fprintf(stderr, "unknown BA method \"%d\" in sba_driver()!\n", howto);
      exit(1);
  }
  end_time=clock();
  

 	fflush(stdout);
 
  fprintf(stdout, "SBA using %d 3D pts, %d frames and %d image projections, %d variables\n", numpts3D, nframes, numprojs, nvars);
  fprintf(stdout, "\nMethod %s, %s driver, %s jacobian\n\n", howtoname[howto],
                  expert? "expert" : "simple",
                  analyticjac? "analytic" : "approximate");
  fprintf(stdout, "SBA returned %d in %g iter, reason %g, error %g [initial %g], %d/%d func/fjac evals, %d lin. systems\n", n, info[5], info[6], info[1]/numprojs, info[0]/numprojs, (int)info[7], (int)info[8], (int)info[9]);
  fprintf(stdout, "Elapsed time: %.2lf seconds, %.2lf msecs\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC,
                  ((double) (end_time - start_time)) / CLOCKS_PER_MSEC);
  fflush(stdout);

  /* refined motion and structure are now in motstruct */

  printSBAData(motstruct, nframes, numpts3D, imgpts, numprojs, vmask); //Min modified 

  /* just in case... */
  globs.intrcalib=NULL;

  free(motstruct);
  free(imgpts);
  free(vmask);
}

int main(int argc, char *argv[])
{
  if(argc!=4){
    fprintf(stderr, "Usage is %s <camera params> <point params> <intrinsic calibration params>\n", argv[0]);
    exit(1);
  }

  sba_driver(argv[1], argv[2], argv[3]);

  return 0;
}
