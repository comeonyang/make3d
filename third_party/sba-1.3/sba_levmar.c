/////////////////////////////////////////////////////////////////////////////////
//// 
////  Expert drivers for sparse bundle adjustment based on the
////  Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "sba.h"
#include "sba_chkjac.h"

#define SBA_EPSILON       1E-12
#define SBA_EPSILON_SQ    ( (SBA_EPSILON)*(SBA_EPSILON) )

#define SBA_ONE_THIRD     0.3333333334 /* 1.0/3.0 */


#define emalloc(sz)       emalloc_(__FILE__, __LINE__, sz)

#define FABS(x)           (((x)>=0)? (x) : -(x))

#define ROW_MAJOR         0
#define COLUMN_MAJOR      1
#define MAT_STORAGE       COLUMN_MAJOR


/* inline */
#ifdef _MSC_VER
#define inline __inline //MSVC
#elif !defined(__GNUC__)
#define inline //other than MSVC, GCC: define empty
#endif

/* contains information necessary for computing a finite difference approximation to a jacobian,
 * e.g. function to differentiate, problem dimensions and pointers to working memory buffers
 */
struct fdj_data_x_ {
  void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata); /* function to differentiate */
  int cnp, pnp, mnp;  /* parameter numbers */
  int *func_rcidxs,
      *func_rcsubs;   /* working memory for func invocations.
                       * Notice that this has to be different
                       * than the working memory used for
                       * evaluating the jacobian!
                       */
  double *hx, *hxx;   /* memory to save results in */
  void *adata;
};

/* auxiliary memory allocation routine with error checking */
inline static void *emalloc_(char *file, int line, size_t sz)
{
void *ptr;

  ptr=(void *)malloc(sz);
  if(ptr==NULL){
    fprintf(stderr, "memory allocation request for %u bytes failed in file %s, line %d, exiting", sz, file, line);
    exit(1);
  }

  return ptr;
}

/* auxiliary routine for clearing an array of doubles */
inline static void _dblzero(register double *arr, register int count)
{
  while(--count>=0)
    *arr++=0.0;
}

/* auxiliary routine for computing the mean reprojection error; used for debugging */
static double sba_mean_repr_error(int n, int mnp, double *x, double *hx, struct sba_crsm *idxij, int *rcidxs, int *rcsubs)
{
register int i, j;
int nnz, nprojs;
double *ptr1, *ptr2;
double err;

  for(i=nprojs=0, err=0.0; i<n; ++i){
    nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
    nprojs+=nnz;
    for(j=0; j<nnz; ++j){ /* point i projecting on camera rcsubs[j] */
      ptr1=x + idxij->val[rcidxs[j]]*mnp;
      ptr2=hx + idxij->val[rcidxs[j]]*mnp;

      err+=sqrt((ptr1[0]-ptr2[0])*(ptr1[0]-ptr2[0]) + (ptr1[1]-ptr2[1])*(ptr1[1]-ptr2[1]));
    }
  }

  return err/((double)(nprojs));
}

/* print the solution in p using sba's text format. If cnp/pnp==0 only points/cameras are printed */
static void sba_print_sol(int n, int m, double *p, int cnp, int pnp, double *x, int mnp, struct sba_crsm *idxij, int *rcidxs, int *rcsubs)
{
register int i, j, ii;
int nnz;
double *ptr;

  if(cnp){
    /* print camera parameters */
    for(j=0; j<m; ++j){
      ptr=p+cnp*j;
      for(ii=0; ii<cnp; ++ii)
        printf("%g ", ptr[ii]);
      printf("\n");
    }
  }

  if(pnp){
    /* 3D & 2D point parameters */
    printf("\n\n\n# X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ...\n");
    for(i=0; i<n; ++i){
      ptr=p+cnp*m+i*pnp;
      for(ii=0; ii<pnp; ++ii) // print 3D coordinates
        printf("%g ", ptr[ii]);

      nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
      printf("%d ", nnz);

      for(j=0; j<nnz; ++j){ /* point i projecting on camera rcsubs[j] */
        ptr=x + idxij->val[rcidxs[j]]*mnp;

        printf("%d ", rcsubs[j]);
        for(ii=0; ii<mnp; ++ii) // print 2D coordinates
          printf("%g ", ptr[ii]);
      }
      printf("\n");
    }
  }
}


/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is approximated with the aid of finite differences and is returned in the order
 * (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/da_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 * Problem-specific information is assumed to be stored in a structure pointed to by "dat".
 *
 * NOTE: The jacobian (for n=4, m=3) in matrix form has the following structure:
 *       A_11  0     0     B_11 0    0    0
 *       0     A_12  0     B_12 0    0    0
 *       0     0     A_13  B_13 0    0    0
 *       A_21  0     0     0    B_21 0    0
 *       0     A_22  0     0    B_22 0    0
 *       0     0     A_23  0    B_23 0    0
 *       A_31  0     0     0    0    B_31 0
 *       0     A_32  0     0    0    B_32 0
 *       0     0     A_33  0    0    B_33 0
 *       A_41  0     0     0    0    0    B_41
 *       0     A_42  0     0    0    0    B_42
 *       0     0     A_43  0    0    0    B_43
 *       To reduce the total number of objective function evaluations, this structure can be
 *       exploited as follows: A certain d is added to the j-th parameters of all cameras and
 *       the objective function is evaluated at the resulting point. This evaluation suffices
 *       to compute the corresponding columns of *all* A_ij through finite differences. A similar
 *       strategy allows the computation of the B_ij. Overall, only cnp+pnp+1 objective function
 *       evaluations are needed to compute the jacobian, much fewer compared to the m*cnp+n*pnp+1
 *       that would be required by the naive strategy of computing one column of the jacobian
 *       per function evaluation. See Nocedal-Wright, ch. 7, pp. 169. Although this approach is
 *       much faster compared to the naive strategy, it is not preferable to analytic jacobians,
 *       since the latter are considerably faster to compute and result in fewer LM iterations.
 */
static void sba_fdjac_x(
    double *p,                /* I: current parameter estimate, (m*cnp+n*pnp)x1 */
    struct sba_crsm *idxij,   /* I: sparse matrix containing the location of x_ij in hx */
    int    *rcidxs,           /* work array for the indexes of nonzero elements of a single sparse matrix row/column */
    int    *rcsubs,           /* work array for the subscripts of nonzero elements in a single sparse matrix row/column */
    double *jac,              /* O: array for storing the approximated jacobian */
    void   *dat)              /* I: points to a "fdj_data_x_" structure */
{
  register int i, j, ii, jj;
  double *pa, *pb, *pqr, *ppt, *jaca, *jacb;
  register double *pAB, *phx, *phxx;
  int n, m, nm, nnz, Asz, Bsz, idx;

  double *tmpd;
  register double d;

  struct fdj_data_x_ *fdjd;
  void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
  double *hx, *hxx;
  int cnp, pnp, mnp;
  void *adata;


  /* retrieve problem-specific information passed in *dat */
  fdjd=(struct fdj_data_x_ *)dat;
  func=fdjd->func;
  cnp=fdjd->cnp; pnp=fdjd->pnp; mnp=fdjd->mnp;
  hx=fdjd->hx;
  hxx=fdjd->hxx;
  adata=fdjd->adata;

  n=idxij->nr; m=idxij->nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp;
  jaca=jac; jacb=jac+idxij->nnz*Asz;

  nm=(n>=m)? n : m; // max(n, m);
  tmpd=(double *)emalloc(nm*sizeof(double));

  (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hx, adata); // evaluate supplied function on current solution

  if(cnp){ // is motion varying?
    /* compute A_ij */
    for(jj=0; jj<cnp; ++jj){
      for(j=0; j<m; ++j){
        pqr=pa+j*cnp; // j-th camera parameters
        /* determine d=max(SBA_DELTA_SCALE*|pqr[jj]|, SBA_MIN_DELTA), see HZ */
        d=(double)(SBA_DELTA_SCALE)*pqr[jj]; // force evaluation
        d=FABS(d);
        if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;

        tmpd[j]=d;
        pqr[jj]+=d;
      }

      (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);

      for(j=0; j<m; ++j){
        pqr=pa+j*cnp; // j-th camera parameters
        pqr[jj]-=tmpd[j]; /* restore */
        d=1.0/tmpd[j]; /* invert so that divisions can be carried out faster as multiplications */

        nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
        for(i=0; i<nnz; ++i){
          idx=idxij->val[rcidxs[i]];
          phx=hx + idx*mnp; // set phx to point to hx_ij
          phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
          pAB=jaca + idx*Asz; // set pAB to point to A_ij

          for(ii=0; ii<mnp; ++ii)
            pAB[ii*cnp+jj]=(phxx[ii]-phx[ii])*d;
        }
      }
    }
  }

  if(pnp){ // is structure varying?
    /* compute B_ij */
    for(jj=0; jj<pnp; ++jj){
      for(i=0; i<n; ++i){
        ppt=pb+i*pnp; // i-th point parameters
        /* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
        d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
        d=FABS(d);
        if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;

        tmpd[i]=d;
        ppt[jj]+=d;
      }

      (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);

      for(i=0; i<n; ++i){
        ppt=pb+i*pnp; // i-th point parameters
        ppt[jj]-=tmpd[i]; /* restore */
        d=1.0/tmpd[i]; /* invert so that divisions can be carried out faster as multiplications */

        nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
        for(j=0; j<nnz; ++j){
          idx=idxij->val[rcidxs[j]];
          phx=hx + idx*mnp; // set phx to point to hx_ij
          phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
          pAB=jacb + idx*Bsz; // set pAB to point to B_ij

          for(ii=0; ii<mnp; ++ii)
            pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
        }
      }
    }
  }

  free(tmpd);
}

/* Bundle adjustment on camera and structure parameters 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 */

int sba_motstr_levmar_x(
    const int n,   /* number of points */
    const int m,   /* number of images */
    const int mcon,/* number of images (starting from the 1st) whose parameters should not be modified.
					          * All A_ij (see below) with j<mcon are assumed to be zero
					          */
    char *vmask,  /* visibility mask: vmask[i][j]=1 if point i visible in image j, 0 otherwise. nxm */
    double *p,    /* initial parameter vector p0: (a1, ..., am, b1, ..., bn).
                   * aj are the image j parameters, bi are the i-th point parameters,
                   * size m*cnp + n*pnp
                   */
    const int cnp,/* number of parameters for ONE camera; e.g. 6 for Euclidean cameras */
    const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
    double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                   * x_ij is the projection of the i-th point on the j-th image.
                   * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                   * see vmask[i][j], max. size n*m*mnp
                   */
    const int mnp,/* number of parameters for EACH measurement; usually 2 */
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                                              /* functional relation describing measurements. Given a parameter vector p,
                                               * computes a prediction of the measurements \hat{x}. p is (m*cnp + n*pnp)x1,
                                               * \hat{x} is (n*m*mnp)x1, maximum
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               */
    void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                                              /* function to evaluate the sparse jacobian dX/dp.
                                               * The Jacobian is returned in jac as
                                               * (dx_11/da_1, ..., dx_1m/da_m, ..., dx_n1/da_1, ..., dx_nm/da_m,
                                               *  dx_11/db_1, ..., dx_1m/db_1, ..., dx_n1/db_n, ..., dx_nm/db_n), or (using HZ's notation),
                                               * jac=(A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm)
                                               * Notice that depending on idxij, some of the A_ij and B_ij might be missing.
                                               * Note also that A_ij and B_ij are mnp x cnp and mnp x pnp matrices resp. and they
                                               * should be stored in jac in row-major order.
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               *
                                               * If NULL, the jacobian is approximated by repetitive func calls and finite
                                               * differences. This is computationally inefficient and thus NOT recommended.
                                               */
    void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

    int itmax,         /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
    int verbose,       /* I: verbosity */
    double opts[SBA_OPTSSZ],
	                     /* I: minim. options [\mu, \epsilon1, \epsilon2]. Respectively the scale factor for initial \mu,
                        * stopping thresholds for ||J^T e||_inf, ||dp||_2 and ||e||_2
                        */
    double info[SBA_INFOSZ]
	                     /* O: information regarding the minimization. Set to NULL if don't care
                        * info[0]=||e||_2 at initial p.
                        * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                        * info[5]= # iterations,
                        * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                        *                                 2 - stopped by small dp
                        *                                 3 - stopped by itmax
                        *                                 4 - singular matrix. Restart from current p with increased mu 
                        *                                 5 - stopped by small ||e||_2
                        *                                 6 - too many attempts to increase damping. Restart with increased mu
                        * info[7]= # function evaluations
                        * info[8]= # jacobian evaluations
			                  * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                        */
)
{
register int i, j, ii, jj, k, l;
int nvis, nnz, retval;

/* The following are work arrays that are dynamically allocated by sba_motstr_levmar_x() */
double *jac;  /* work array for storing the jacobian, max. size n*m*(mnp*cnp + mnp*pnp) */
double *U;    /* work array for storing the U_j in the order U_1, ..., U_m, size m*cnp*cnp */
double *V;    /* work array for storing the V_i in the order V_1, ..., V_n, size n*pnp*pnp */
double *V_1;  /* work array for storing the (V*_i)^-1 in the order (V*_1)^-1, ..., (V*_n)^-1, size n*pnp*pnp */

double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                 max. size n*m*mnp */
double *eab;  /* work array for storing the ea_j & eb_i in the order ea_1, .. ea_m eb_1, .. eb_n size m*cnp + n*pnp */

double *E;   /* work array for storing the e_j in the order e_1, .. e_m, size m*cnp */

/* Notice that the blocks W_ij, Y_ij are zero iff A_ij (equivalently B_ij) is zero. This means
 * that the matrices consisting of blocks W_ij and Y_ij are themselves sparse, similarly to the
 * block matrices made up of the A_ij and B_ij (i.e. jaca and jacb)
 */
double *W;    /* work array for storing the W_ij in the order W_11, ..., W_1m, ..., W_n1, ..., W_nm,
                 max. size n*m*cnp*pnp */
double *Y;    /* work array for storing the Y_ij in the order Y_11, ..., Y_1m, ..., Y_n1, ..., Y_nm,
                 max. size n*m*cnp*pnp */
double *YWt;  /* work array for storing \sum_i Y_ij W_ik^T, size cnp*cnp */
double *S;    /* work array for storing the block array S_jk, size m*m*cnp*cnp */
double *dp;   /* work array for storing the parameter vector updates da_1, ..., da_m, db_1, ..., db_n, size m*cnp + n*pnp */
double *Wtda; /* work array for storing \sum_j W_ij^T da_j, size pnp */

/* Of the above arrays, jac, e, W, Y are sparse and
 * U, V, V_1, eab, E, S, dp are dense. Sparse arrays are indexed through
 * idxij (see below), that is with the same mechanism as the input 
 * measurements vector x
 */

double *pa, *pb, *jaca, *jacb, *ea, *eb, *dpa, *dpb; /* pointers into p, jac, eab and dp respectively */

/* submatrices sizes */
int Asz, Bsz, Usz, Vsz,
    Wsz, Ysz, esz, easz, ebsz,
    YWtsz, Wtdasz, Sblsz;

int Sdim; /* S matrix actual dimension */

register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location of A_ij 
                        * in jaca, B_ij in jacb, etc.
                        * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                        * index k and it is used to efficiently lookup the memory locations where the non-zero
                        * blocks of a sparse matrix/vector are stored
                        */
int maxnm=(n>=m)? n:m, /* max. of (n, m) */
    *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                 column in a sparse matrix, size max(n, m) */
    *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                 sparse matrix, size max(n, m) */

/* The following variables are needed by the LM algorithm */
register int itno;  /* iteration counter */
int issolved;
/* temporary work arrays that are dynamically allocated */
double *hx,         /* \hat{x}_i, max. size m*n*mnp */
       *diagUV,     /* diagonals of U_j, V_i, size m*cnp + n*pnp */
       *pdp;        /* p + dp, size m*cnp + n*pnp */

double *diagU, *diagV; /* pointers into diagUV */

register double mu,  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
double p_eL2, eab_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
double p_L2, dp_L2=DBL_MAX, dF, dL;
double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2], eps3_sq=opts[3]*opts[3];
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0, nlss=0;
int nobs, nvars;
const int mmcon=m-mcon;

struct fdj_data_x_ fdj_data;
void *jac_adata;

/* Initialization */

  /* block sizes */
  Asz=mnp * cnp; Bsz=mnp * pnp;
  Usz=cnp * cnp; Vsz=pnp * pnp;
  Wsz=cnp * pnp; Ysz=cnp * pnp;
  esz=mnp;
  easz=cnp; ebsz=pnp;
  YWtsz=cnp * cnp;
  Wtdasz=pnp;
  Sblsz=cnp * cnp;
  Sdim=mmcon * cnp;

  /* count total number of visible image points */
  for(i=nvis=0, jj=n*m; i<jj; ++i)
    nvis+=vmask[i];

  nobs=nvis*mnp;
  nvars=m*cnp + n*pnp;
  if(nobs<nvars){
    fprintf(stderr, "sba_motstr_levmar_x(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
    exit(1);
  }

  /* allocate work arrays */
  jac=(double *)emalloc(nvis*(Asz+Bsz)*sizeof(double));
  U=(double *)emalloc(m*Usz*sizeof(double));
  V=(double *)emalloc(n*Vsz*sizeof(double));
  V_1=(double *)emalloc(n*Vsz*sizeof(double));
  e=(double *)emalloc(nobs*sizeof(double));
  eab=(double *)emalloc(nvars*sizeof(double));
  E=(double *)emalloc(m*cnp*sizeof(double));
  W=(double *)emalloc(nvis*Wsz*sizeof(double));
  Y=(double *)emalloc(nvis*Ysz*sizeof(double));
  YWt=(double *)emalloc(YWtsz*sizeof(double));
  S=(double *)emalloc(m*m*Sblsz*sizeof(double));
  dp=(double *)emalloc(nvars*sizeof(double));
  Wtda=(double *)emalloc(pnp*sizeof(double));
  rcidxs=(int *)emalloc(maxnm*sizeof(int));
  rcsubs=(int *)emalloc(maxnm*sizeof(int));


  hx=(double *)emalloc(nobs*sizeof(double));
  diagUV=(double *)emalloc(nvars*sizeof(double));
  pdp=(double *)emalloc(nvars*sizeof(double));


  /* set up auxiliary pointers */
  pa=p; pb=p+m*cnp;
  jaca=jac; jacb=jac+nvis*Asz;
  ea=eab; eb=eab+m*cnp;
  dpa=dp; dpb=dp+m*cnp;

  diagU=diagUV; diagV=diagUV + m*cnp;


  /* allocate & fill up the idxij structure */
  sba_crsm_alloc(&idxij, n, m, nvis);
  for(i=k=0; i<n; ++i){
    idxij.rowptr[i]=k;
    ii=i*m;
    for(j=0; j<m; ++j)
      if(vmask[ii+j]){
        idxij.val[k]=k;
        idxij.colidx[k++]=j;
      }
  }
  idxij.rowptr[n]=nvis;

  /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
  if(!fjac){
    fdj_data.func=func;
    fdj_data.cnp=cnp;
    fdj_data.pnp=pnp;
    fdj_data.mnp=mnp;
    fdj_data.hx=hx;
    fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
    fdj_data.func_rcidxs=(int *)emalloc(2*maxnm*sizeof(int));
    fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxnm;
    fdj_data.adata=adata;

    fjac=sba_fdjac_x;
    jac_adata=(void *)&fdj_data;
  }
  else{
    fdj_data.hxx=NULL;
    jac_adata=adata;
  }

  if(itmax==0){ /* verify jacobian */
    sba_motstr_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, mcon, cnp, pnp, mnp, adata, jac_adata);
    retval=0;
    goto freemem_and_return;
  }

  /* compute the error vectors e_ij in hx */
  (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
  /* compute e=x - f(p) and its L2 norm */
  for(i=0, p_eL2=0.0; i<nobs; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }

  if(verbose) printf("initial motstr-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
  init_p_eL2=p_eL2;

  for(itno=stop=0; itno<itmax && !stop; ++itno){
    /* Note that p, e and ||e||_2 have been updated at the previous iteration */

    /* compute derivative submatrices A_ij, B_ij */
    (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

    /* compute U_j = \sum_i A_ij^T A_ij */ // \Sigma here!
    /* U_j is symmetric, therefore its computation can be speeded up by
     * computing only the upper part and then reusing it for the lower one.
     * Recall that A_ij is mnp x cnp
     */
    /* Also compute ea_j = \sum_i A_ij^T e_ij */ // \Sigma here!
    /* Recall that e_ij is mnp x 1
     */
    _dblzero(U, m*Usz); /* clear all U_j */
    _dblzero(ea, m*easz); /* clear all ea_j */
    for(j=mcon; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      ptr2=ea + j*easz; // set ptr2 to point to ea_j

      nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
      for(i=0; i<nnz; ++i){
        /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
        ptr3=jaca + idxij.val[rcidxs[i]]*Asz;

        /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
        for(ii=0; ii<cnp; ++ii){
          for(jj=ii; jj<cnp; ++jj){
            for(k=0, sum=0.0; k<mnp; ++k)
              sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
            ptr1[ii*cnp+jj]+=sum;
          }

          /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
          for(jj=0; jj<ii; ++jj)
            ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
        }

        ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
        /* compute A_ij^T e_ij and add it to ea_j */
        for(ii=0; ii<cnp; ++ii){
          for(jj=0, sum=0.0; jj<mnp; ++jj)
            sum+=ptr3[jj*cnp+ii]*ptr4[jj];
          ptr2[ii]+=sum;
        }
      }
    }

    /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
    /* V_i is symmetric, therefore its computation can be speeded up by
     * computing only the upper part and then reusing it for the lower one.
     * Recall that B_ij is mnp x pnp
     */
    /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
    /* Recall that e_ij is mnp x 1
     */
	  _dblzero(V, n*Vsz); /* clear all V_i */
	  _dblzero(eb, n*ebsz); /* clear all eb_i */
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      ptr2=eb + i*ebsz; // set ptr2 to point to eb_i

      nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
      for(j=0; j<nnz; ++j){
        /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
        ptr3=jacb + idxij.val[rcidxs[j]]*Bsz;
      
        /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
        for(ii=0; ii<pnp; ++ii){
          for(jj=ii; jj<pnp; ++jj){
            for(k=0, sum=0.0; k<mnp; ++k)
              sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
            ptr1[ii*pnp+jj]+=sum;
          }

          /* copy the LOWER TRIANGULAR PART of V_i from the upper one */
          for(jj=0; jj<ii; ++jj)
            ptr1[ii*pnp+jj]=ptr1[jj*pnp+ii];
        }

        ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
        /* compute B_ij^T e_ij and add it to eb_i */
        for(ii=0; ii<pnp; ++ii){
          for(jj=0, sum=0.0; jj<mnp; ++jj)
            sum+=ptr3[jj*pnp+ii]*ptr4[jj];
          ptr2[ii]+=sum;
        }
      }
    }

    /* compute W_ij =  A_ij^T B_ij */ // \Sigma here!
    /* Recall that A_ij is mnp x cnp and B_ij is mnp x pnp
     */
    _dblzero(W, nvis*Wsz); /* clear all W_ij */
    for(i=0; i<n; ++i){
      nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
      for(j=0; j<nnz; ++j){
        /* set ptr1 to point to W_ij, actual column number in rcsubs[j] */
        if(rcsubs[j]<mcon) continue; /* A_ij is zero */

        ptr1=W + idxij.val[rcidxs[j]]*Wsz;
        /* set ptr2 & ptr3 to point to A_ij & B_ij resp. */
        ptr2=jaca + idxij.val[rcidxs[j]]*Asz;
        ptr3=jacb + idxij.val[rcidxs[j]]*Bsz;
        /* compute A_ij^T B_ij and store it in W_ij */
        for(ii=0; ii<cnp; ++ii)
          for(jj=0; jj<pnp; ++jj){
            for(k=0, sum=0.0; k<mnp; ++k)
              sum+=ptr2[k*cnp+ii]*ptr3[k*pnp+jj];
            ptr1[ii*pnp+jj]=sum;
          }
      }
    }

    /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=0, p_L2=eab_inf=0.0; i<nvars; ++i){
      if(eab_inf < (tmp=FABS(eab[i]))) eab_inf=tmp;
      p_L2+=p[i]*p[i];
    }
    //p_L2=sqrt(p_L2);

    /* save diagonal entries so that augmentation can be later canceled.
     * Diagonal entries are in U_j and V_i
     */
    for(j=mcon; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
      for(i=0; i<cnp; ++i)
        ptr2[i]=ptr1[i*cnp+i];
    }
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
      for(j=0; j<pnp; ++j)
        ptr2[j]=ptr1[j*pnp+j];
    }

/*
if(!(itno%100)){
  printf("Current estimate: ");
  for(i=0; i<nvars; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g\n", eab_inf, p_eL2);
}
*/

    /* check for convergence */
    if((eab_inf <= eps1)){
      dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(itno==0){
      for(i=mcon*cnp, tmp=DBL_MIN; i<nvars; ++i)
        if(diagUV[i]>tmp) tmp=diagUV[i]; /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */
    while(1){
      /* augment U, V */
      for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        for(i=0; i<cnp; ++i)
          ptr1[i*cnp+i]+=mu;
      }
      for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        for(j=0; j<pnp; ++j)
          ptr1[j*pnp+j]+=mu;

		    /* compute V*_i^-1 */
		    ptr2=V_1 + i*Vsz; // set ptr2 to point to (V*_i)^-1
        /* inverting V*_i with LU seems to be faster than Cholesky */
        j=sba_mat_invert_LU(ptr1, ptr2, pnp); //invert the pnp x pnp matrix pointed to by ptr1 and save it in ptr2
        //j=sba_mat_invert_Chol(ptr1, ptr2, pnp); //invert the pnp x pnp matrix pointed to by ptr1 and save it in ptr2
		    if(!j){
			    fprintf(stderr, "Singular matrix V*_i (i=%d) in sba_motstr_levmar_x()\n", i);
			    exit(1);
		    }
      }

      /* compute Y_ij = W_ij (V*_i)^-1
       * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
       */
      _dblzero(Y, nvis*Ysz); /* clear all Y_ij */
      for(i=0; i<n; ++i){
        /* set ptr3 to point to (V*_i)^-1 */
        ptr3=V_1 + i*Vsz;
        nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
        for(j=0; j<nnz; ++j){
          /* set ptr1 to point to Y_ij, actual column number in rcsubs[j] */
		      if(rcsubs[j]<mcon) continue; /* W_ij is zero */

          ptr1=Y + idxij.val[rcidxs[j]]*Ysz;
          /* set ptr2 to point to W_ij resp. */
          ptr2=W + idxij.val[rcidxs[j]]*Wsz;
          /* compute W_ij (V*_i)^-1 and store it in Y_ij */
          for(ii=0; ii<cnp; ++ii)
            for(jj=0; jj<pnp; ++jj){
              for(k=0, sum=0.0; k<pnp; ++k)
                sum+=ptr2[ii*pnp+k]*ptr3[k*pnp+jj];
              ptr1[ii*pnp+jj]=sum;
            }
        }
      }

      _dblzero(E, m*easz); /* clear all e_j */
      /* compute the mmcon x mmcon block matrix S and e_j */

      /* Note that S is symmetric, therefore its computation can be
       * speeded up by computing only the upper part and then reusing
       * it for the lower one.
		   */
      for(j=mcon; j<m; ++j){
		    nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero Y_ij, i=0...n-1 */

        /* compute the UPPER TRIANGULAR PART of S */
        for(k=j; k<m; ++k){ // j>=mcon
          /* compute \sum_i Y_ij W_ik^T in YWt */
          /* Recall that Y_ij is cnp x pnp and W_ik is cnp x pnp */ 
          _dblzero(YWt, YWtsz); /* clear YWt */

          for(i=0; i<nnz; ++i){
            register double *pYWt;

            /* find the min and max column indices of the elements in row i (actually rcsubs[i])
             * and make sure that k falls within them. This test handles W_ik's which are
             * certain to be zero without bothering to call sba_crsm_elmidx()
             */
            ii=idxij.colidx[idxij.rowptr[rcsubs[i]]];
            jj=idxij.colidx[idxij.rowptr[rcsubs[i]+1]-1];
            if(k<ii || k>jj) continue; /* W_ik == 0 */

            /* set ptr2 to point to W_ik */
            l=sba_crsm_elmidx(&idxij, rcsubs[i], k);
            if(l==-1) continue; /* W_ik == 0 */

            ptr2=W + idxij.val[l]*Wsz;
            /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
            ptr1=Y + idxij.val[rcidxs[i]]*Ysz;
            for(ii=0; ii<cnp; ++ii){
              ptr3=ptr1+ii*pnp;
              pYWt=YWt+ii*cnp;

              for(jj=0; jj<cnp; ++jj){
                ptr4=ptr2+jj*pnp;
                for(l=0, sum=0.0; l<pnp; ++l)
                  sum+=ptr3[l]*ptr4[l]; //ptr1[ii*pnp+l]*ptr2[jj*pnp+l];
                pYWt[jj]+=sum; //YWt[ii*cnp+jj]+=sum;
              }
            }
          }

          ptr1=U + j*Usz; // set ptr1 to point to U_j
		  
		      /* since the linear system involving S is solved with lapack,
		       * it is preferable to store S in column major (i.e. fortran)
		       * order, so as to avoid unecessary transposing/copying.
           */
#if MAT_STORAGE==COLUMN_MAJOR
          ptr2=S + (k-mcon)*mmcon*Usz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#else
          ptr2=S + (j-mcon)*mmcon*Usz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#endif
		  
          if(j!=k){ /* Kronecker */
            for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
              for(jj=0; jj<cnp; ++jj)
                ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
				                -YWt[jj*cnp+ii];
#else
				                -YWt[ii*cnp+jj];
#endif
          }
          else{
            for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
              for(jj=0; jj<cnp; ++jj)
                ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
				                ptr1[jj*cnp+ii] - YWt[jj*cnp+ii];
#else
				                ptr1[ii*cnp+jj] - YWt[ii*cnp+jj];
#endif
          }
        }

        /* copy the LOWER TRIANGULAR PART of S from the upper one */
        for(k=mcon; k<j; ++k){
#if MAT_STORAGE==COLUMN_MAJOR
          ptr1=S + (k-mcon)*mmcon*Usz + (j-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
          ptr2=S + (j-mcon)*mmcon*Usz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#else
          ptr1=S + (j-mcon)*mmcon*Usz + (k-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
          ptr2=S + (k-mcon)*mmcon*Usz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#endif
          for(ii=0; ii<cnp; ++ii, ptr1+=Sdim)
            for(jj=0, ptr3=ptr2+ii; jj<cnp; ++jj, ptr3+=Sdim)
              ptr1[jj]=*ptr3;
        }

        /* compute e_j=ea_j - \sum_i Y_ij eb_i */
        /* Recall that Y_ij is cnp x pnp and eb_i is pnp x 1 */
        ptr1=E + j*easz; // set ptr1 to point to e_j

        for(i=0; i<nnz; ++i){
          /* set ptr2 to point to Y_ij, actual row number in rcsubs[i] */
          ptr2=Y + idxij.val[rcidxs[i]]*Ysz;

          /* set ptr3 to point to eb_i */
          ptr3=eb + rcsubs[i]*ebsz;
          for(ii=0; ii<cnp; ++ii){
            ptr4=ptr2+ii*pnp;
            for(jj=0, sum=0.0; jj<pnp; ++jj)
              sum+=ptr4[jj]*ptr3[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
            ptr1[ii]+=sum;
          }
        }

        ptr2=ea + j*easz; // set ptr2 to point to ea_j
        for(i=0; i<easz; ++i)
          ptr1[i]=ptr2[i] - ptr1[i];
      }

#if 0
      if(verbose>1){ /* count the nonzeros in S */
        for(i=ii=0; i<Sdim*Sdim; ++i)
          if(S[i]!=0.0) ++ii;
        printf("\nS density %10g\n", ((double)ii)/(Sdim*Sdim));

      }
#endif

      /* solve the linear system S dpa = E to compute the da_j.
       *
       * Note that if MAT_STORAGE==1 S is modified in the following call;
       * this is OK since S is recomputed for each iteration
       */
	    //issolved=sba_Axb_LU(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); ++nlss;
      issolved=sba_Axb_Chol(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); ++nlss;
      //issolved=sba_Axb_QRnoQ(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); ++nlss;
      //issolved=sba_Axb_QR(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); ++nlss;
	    //issolved=sba_Axb_SVD(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); ++nlss;
	    //issolved=sba_Axb_CG(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, (3*Sdim)/2, 1E-10, SBA_CG_JACOBI, MAT_STORAGE); ++nlss;

	    _dblzero(dpa, mcon*cnp); /* no change for the first mcon camera params */

      if(issolved){

        /* compute the db_i */
        for(i=0; i<n; ++i){
          ptr1=dpb + i*ebsz; // set ptr1 to point to db_i

          /* compute \sum_j W_ij^T da_j */
          /* Recall that W_ij is cnp x pnp and da_j is cnp x 1 */
          _dblzero(Wtda, Wtdasz); /* clear Wtda */
          nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
          for(j=0; j<nnz; ++j){
            /* set ptr2 to point to W_ij, actual column number in rcsubs[j] */
			      if(rcsubs[j]<mcon) continue; /* W_ij is zero */

            ptr2=W + idxij.val[rcidxs[j]]*Wsz;

            /* set ptr3 to point to da_j */
            ptr3=dpa + rcsubs[j]*cnp;

            for(ii=0; ii<pnp; ++ii){
              ptr4=ptr2+ii;
              for(jj=0, sum=0.0; jj<cnp; ++jj)
                sum+=ptr4[jj*pnp]*ptr3[jj]; //ptr2[jj*pnp+ii]*ptr3[jj];
              Wtda[ii]+=sum;
            }
          }

          /* compute eb_i - \sum_j W_ij^T da_j = eb_i - Wtda in Wtda */
          ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
          for(ii=0; ii<pnp; ++ii)
            Wtda[ii]=ptr2[ii] - Wtda[ii];

          /* compute the product (V*_i)^-1 Wtda = (V*_i)^-1 (eb_i - \sum_j W_ij^T da_j) */
          ptr2=V_1 + i*Vsz; // set ptr2 to point to (V*_i)^-1
          for(ii=0; ii<pnp; ++ii){
            for(jj=0, sum=0.0; jj<pnp; ++jj)
              sum+=ptr2[ii*pnp+jj]*Wtda[jj];
            ptr1[ii]=sum;
          }
        }

        /* parameter vector updates are now in dpa, dpb */

        /* compute p's new estimate and ||dp||^2 */
        for(i=0, dp_L2=0.0; i<nvars; ++i){
          pdp[i]=p[i] + (tmp=dp[i]);
          dp_L2+=tmp*tmp;
        }
        //dp_L2=sqrt(dp_L2);

        if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
       //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
         stop=4;
         break;
       }

        (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
        if(verbose>1)
          printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
        for(i=0, pdp_eL2=0.0; i<nobs; ++i){ /* compute ||e(pdp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pdp_eL2+=tmp*tmp;
        }

        for(i=0, dL=0.0; i<nvars; ++i)
          dL+=dp[i]*(mu*dp[i]+eab[i]);

        dF=p_eL2-pdp_eL2;

        if(verbose>1)
          printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(2.0*dF/dL-1.0);
          tmp=1.0-tmp*tmp*tmp;
          mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
          nu=2;

          for(i=0; i<nvars; ++i) /* update p's estimate */
            p[i]=pdp[i];

          for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pdp_eL2;
          break;
        }
      } /* issolved */

      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown) */
        fprintf(stderr, "Too many failed attempts to increase the damping factor in sba_motstr_levmar_x()! Singular hessian matrix?\n");
        //exit(1);
        stop=6;
        break;
      }
      nu=nu2;

      /* restore U, V diagonal entries */
      for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
        for(i=0; i<cnp; ++i)
          ptr1[i*cnp+i]=ptr2[i];
      }
      for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
        for(j=0; j<pnp; ++j)
          ptr1[j*pnp+j]=ptr2[j];
      }
    } /* inner loop */

    if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
  }

  if(itno>=itmax) stop=3;

  /* restore U, V diagonal entries */
  for(j=mcon; j<m; ++j){
    ptr1=U + j*Usz; // set ptr1 to point to U_j
    ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
    for(i=0; i<cnp; ++i)
      ptr1[i*cnp+i]=ptr2[i];
  }
  for(i=0; i<n; ++i){
    ptr1=V + i*Vsz; // set ptr1 to point to V_i
    ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
    for(j=0; j<pnp; ++j)
     ptr1[j*pnp+j]=ptr2[j];
  }

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=eab_inf;
    info[3]=dp_L2;
    for(j=mcon, tmp=DBL_MIN; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      for(i=0; i<cnp; ++i)
        if(tmp<ptr1[i*cnp+i]) tmp=ptr1[i*cnp+i];
    }
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      for(j=0; j<pnp; ++j)
        if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
      }
    info[4]=mu/tmp;
    info[5]=itno;
    info[6]=stop;
    info[7]=nfev;
    info[8]=njev;
    info[9]=nlss;
  }
                                                               
  //sba_print_sol(n, m, p, cnp, pnp, x, mnp, &idxij, rcidxs, rcsubs);
  retval=(stop!=4)?  itno : -1;

freemem_and_return: /* NOTE: this point is also reached via a goto! */

   /* free whatever was allocated */
  free(jac); free(U); free(V);
  free(V_1); free(e); free(eab);  
  free(E);   free(W); free(Y);          
  free(YWt); free(S); free(dp);               
  free(Wtda);
  free(rcidxs); free(rcsubs);

  free(hx); free(diagUV); free(pdp);
  if(fdj_data.hxx){ // cleanup
    free(fdj_data.hxx);
    free(fdj_data.func_rcidxs);
  }

  sba_crsm_free(&idxij);

  return retval;
}


/* Bundle adjustment on camera parameters only 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 */

int sba_mot_levmar_x(
    const int n,   /* number of points */
    const int m,   /* number of images */
    const int mcon,/* number of images (starting from the 1st) whose parameters should not be modified.
					          * All A_ij (see below) with j<mcon are assumed to be zero
					          */
    char *vmask,  /* visibility mask: vmask[i][j]=1 if point i visible in image j, 0 otherwise. nxm */
    double *p,    /* initial parameter vector p0: (a1, ..., am).
                   * aj are the image j parameters, size m*cnp */
    const int cnp,/* number of parameters for ONE camera; e.g. 6 for Euclidean cameras */
    double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                   * x_ij is the projection of the i-th point on the j-th image.
                   * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                   * see vmask[i][j], max. size n*m*mnp
                   */
    const int mnp,/* number of parameters for EACH measurement; usually 2 */
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                                              /* functional relation describing measurements. Given a parameter vector p,
                                               * computes a prediction of the measurements \hat{x}. p is (m*cnp)x1,
                                               * \hat{x} is (n*m*mnp)x1, maximum
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               */
    void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                                              /* function to evaluate the sparse jacobian dX/dp.
                                               * The Jacobian is returned in jac as
                                               * (dx_11/da_1, ..., dx_1m/da_m, ..., dx_n1/da_1, ..., dx_nm/da_m), or (using HZ's notation),
                                               * jac=(A_11, ..., A_1m, ..., A_n1, ..., A_nm)
                                               * Notice that depending on idxij, some of the A_ij might be missing.
                                               * Note also that the A_ij are mnp x cnp matrices and they
                                               * should be stored in jac in row-major order.
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               *
                                               * If NULL, the jacobian is approximated by repetitive func calls and finite
                                               * differences. This is computationally inefficient and thus NOT recommended.
                                               */
    void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

    int itmax,         /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
    int verbose,       /* I: verbosity */
    double opts[SBA_OPTSSZ],
	                     /* I: minim. options [\mu, \epsilon1, \epsilon2]. Respectively the scale factor for initial \mu,
                        * stopping thresholds for ||J^T e||_inf, ||dp||_2 and ||e||_2
                        */
    double info[SBA_INFOSZ]
	                     /* O: information regarding the minimization. Set to NULL if don't care
                        * info[0]=||e||_2 at initial p.
                        * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                        * info[5]= # iterations,
                        * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                        *                                 2 - stopped by small dp
                        *                                 3 - stopped by itmax
                        *                                 4 - singular matrix. Restart from current p with increased mu 
                        *                                 5 - stopped by small ||e||_2
                        *                                 6 - too many attempts to increase damping. Restart with increased mu
                        * info[7]= # function evaluations
                        * info[8]= # jacobian evaluations
			                  * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                        */
)
{
register int i, j, ii, jj, k;
int nvis, nnz, retval;

/* The following are work arrays that are dynamically allocated by sba_mot_levmar_x() */
double *jac; /* work array for storing the jacobian, max. size n*m*mnp*cnp */
double *U;    /* work array for storing the U_j in the order U_1, ..., U_m, size m*cnp*cnp */

double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                 max. size n*m*mnp */
double *ea;   /* work array for storing the ea_j in the order ea_1, .. ea_m, size m*cnp */

double *dp;   /* work array for storing the parameter vector updates da_1, ..., da_m, size m*cnp */

/* Of the above arrays, jac, e are sparse and
 * U, ea, dp are dense. Sparse arrays are indexed through
 * idxij (see below), that is with the same mechanism as the input 
 * measurements vector x
 */

/* submatrices sizes */
int Asz, Usz,
    esz, easz;

register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location of A_ij 
                        * in jac, e_ij in e, etc.
                        * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                        * index k and it is used to efficiently lookup the memory locations where the non-zero
                        * blocks of a sparse matrix/vector are stored
                        */
int maxnm=(n>=m)? n:m, /* max. of (n, m) */
    *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                 column in a sparse matrix, size max(n, m) */
    *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                 sparse matrix, size max(n, m) */

/* The following variables are needed by the LM algorithm */
register int itno;  /* iteration counter */
int nsolved;
/* temporary work arrays that are dynamically allocated */
double *hx,         /* \hat{x}_i, max. size m*n*mnp */
       *diagU,      /* diagonals of U_j, size m*cnp */
       *pdp;        /* p + dp, size m*cnp */

register double mu,  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
double p_eL2, ea_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
double p_L2, dp_L2=DBL_MAX, dF, dL;
double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2], eps3_sq=opts[3]*opts[3];
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0, nlss=0;
int nobs, nvars;

struct fdj_data_x_ fdj_data;
void *jac_adata;

/* Initialization */

  /* block sizes */
  Asz=mnp * cnp; Usz=cnp * cnp;
  esz=mnp; easz=cnp;
  
  /* count total number of visible image points */
  for(i=nvis=0, jj=n*m; i<jj; ++i)
    nvis+=vmask[i];

  nobs=nvis*mnp;
  nvars=m*cnp;
  if(nobs<nvars){
    fprintf(stderr, "sba_mot_levmar_x(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
    exit(1);
  }

  /* allocate work arrays */
  jac=(double *)emalloc(nvis*Asz*sizeof(double));
  U=(double *)emalloc(m*Usz*sizeof(double));
  e=(double *)emalloc(nobs*sizeof(double));
  ea=(double *)emalloc(nvars*sizeof(double));
  dp=(double *)emalloc(nvars*sizeof(double));
  rcidxs=(int *)emalloc(maxnm*sizeof(int));
  rcsubs=(int *)emalloc(maxnm*sizeof(int));


  hx=(double *)emalloc(nobs*sizeof(double));
  diagU=(double *)emalloc(nvars*sizeof(double));
  pdp=(double *)emalloc(nvars*sizeof(double));

  /* allocate & fill up the idxij structure */
  sba_crsm_alloc(&idxij, n, m, nvis);
  for(i=k=0; i<n; ++i){
    idxij.rowptr[i]=k;
    ii=i*m;
    for(j=0; j<m; ++j)
      if(vmask[ii+j]){
        idxij.val[k]=k;
        idxij.colidx[k++]=j;
      }
  }
  idxij.rowptr[n]=nvis;

  /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
  if(!fjac){
    fdj_data.func=func;
    fdj_data.cnp=cnp;
    fdj_data.pnp=0;
    fdj_data.mnp=mnp;
    fdj_data.hx=hx;
    fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
    fdj_data.func_rcidxs=(int *)emalloc(2*maxnm*sizeof(int));
    fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxnm;
    fdj_data.adata=adata;

    fjac=sba_fdjac_x;
    jac_adata=(void *)&fdj_data;
  }
  else{
    fdj_data.hxx=NULL;
    jac_adata=adata;
  }

  if(itmax==0){ /* verify jacobian */
    sba_mot_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, mcon, cnp, mnp, adata, jac_adata);
    retval=0;
    goto freemem_and_return;
  }

  /* compute the error vectors e_ij in hx */
  (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
  /* compute e=x - f(p) and its L2 norm */
  for(i=0, p_eL2=0.0; i<nobs; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }

  if(verbose) printf("initial mot-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
  init_p_eL2=p_eL2;

  for(itno=stop=0; itno<itmax && !stop; ++itno){
    /* Note that p, e and ||e||_2 have been updated at the previous iteration */

    /* compute derivative submatrices A_ij */
    (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

    /* compute U_j = \sum_i A_ij^T A_ij */ // \Sigma here!
    /* U_j is symmetric, therefore its computation can be speeded up by
     * computing only the upper part and then reusing it for the lower one.
     * Recall that A_ij is mnp x cnp
     */
    /* Also compute ea_j = \sum_i A_ij^T e_ij */ // \Sigma here!
    /* Recall that e_ij is mnp x 1
     */
	  _dblzero(U, m*Usz); /* clear all U_j */
	  _dblzero(ea, m*easz); /* clear all ea_j */
    for(j=mcon; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      ptr2=ea + j*easz; // set ptr2 to point to ea_j

      nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
      for(i=0; i<nnz; ++i){
        /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
        ptr3=jac + idxij.val[rcidxs[i]]*Asz;

        /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
        for(ii=0; ii<cnp; ++ii){
          for(jj=ii; jj<cnp; ++jj){
            for(k=0, sum=0.0; k<mnp; ++k)
              sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
            ptr1[ii*cnp+jj]+=sum;
          }

          /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
          for(jj=0; jj<ii; ++jj)
            ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
        }

        ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
        /* compute A_ij^T e_ij and add it to ea_j */
        for(ii=0; ii<cnp; ++ii){
          for(jj=0, sum=0.0; jj<mnp; ++jj)
            sum+=ptr3[jj*cnp+ii]*ptr4[jj];
          ptr2[ii]+=sum;
        }
      }
    }

    /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=0, p_L2=ea_inf=0.0; i<nvars; ++i){
      if(ea_inf < (tmp=FABS(ea[i]))) ea_inf=tmp;
      p_L2+=p[i]*p[i];
    }
    //p_L2=sqrt(p_L2);

    /* save diagonal entries so that augmentation can be later canceled.
     * Diagonal entries are in U_j
     */
    for(j=mcon; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
      for(i=0; i<cnp; ++i)
        ptr2[i]=ptr1[i*cnp+i];
    }

/*
if(!(itno%100)){
  printf("Current estimate: ");
  for(i=0; i<nvars; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g\n", ea_inf, p_eL2);
}
*/

    /* check for convergence */
    if((ea_inf <= eps1)){
      dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(itno==0){
      for(i=mcon*cnp, tmp=DBL_MIN; i<nvars; ++i)
        if(diagU[i]>tmp) tmp=diagU[i]; /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */
    while(1){
      /* augment U */
      for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        for(i=0; i<cnp; ++i)
          ptr1[i*cnp+i]+=mu;
      }
 
      /* solve the linear systems U_j da_j = ea_j to compute the da_j */
      _dblzero(dp, mcon*cnp); /* no change for the first mcon camera params */
      for(j=nsolved=mcon; j<m; ++j){
		    ptr1=U + j*Usz; // set ptr1 to point to U_j
		    ptr2=ea + j*easz; // set ptr2 to point to ea_j
		    ptr3=dp + j*cnp; // set ptr3 to point to da_j

		    //nsolved+=sba_Axb_LU(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
		    nsolved+=sba_Axb_Chol(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
		    //nsolved+=sba_Axb_BK(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
		    //nsolved+=sba_Axb_QRnoQ(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
		    //nsolved+=sba_Axb_QR(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
		    //nsolved+=sba_Axb_SVD(ptr1, ptr2, ptr3, cnp, 0); ++nlss;
	      //nsolved+=(sba_Axb_CG(ptr1, ptr2, ptr3, cnp, (3*cnp)/2, 1E-10, SBA_CG_JACOBI, 0) > 0); ++nlss;
	    }

	    if(nsolved==m){

        /* parameter vector updates are now in dp */

        /* compute p's new estimate and ||dp||^2 */
        for(i=0, dp_L2=0.0; i<nvars; ++i){
          pdp[i]=p[i] + (tmp=dp[i]);
          dp_L2+=tmp*tmp;
        }
        //dp_L2=sqrt(dp_L2);

        if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
       //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
         stop=4;
         break;
       }

        (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
        if(verbose>1)
          printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
        for(i=0, pdp_eL2=0.0; i<nobs; ++i){ /* compute ||e(pdp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pdp_eL2+=tmp*tmp;
        }

        for(i=0, dL=0.0; i<nvars; ++i)
          dL+=dp[i]*(mu*dp[i]+ea[i]);

        dF=p_eL2-pdp_eL2;

        if(verbose>1)
          printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(2.0*dF/dL-1.0);
          tmp=1.0-tmp*tmp*tmp;
          mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
          nu=2;

          for(i=0; i<nvars; ++i) /* update p's estimate */
            p[i]=pdp[i];

          for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pdp_eL2;
          break;
        }
      } /* nsolved==m */

      /* if this point is reached, either at least one linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown) */
        fprintf(stderr, "Too many failed attempts to increase the damping factor in sba_mot_levmar_x()! Singular hessian matrix?\n");
        //exit(1);
        stop=6;
        break;
      }
      nu=nu2;

      /* restore U diagonal entries */
      for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
        for(i=0; i<cnp; ++i)
          ptr1[i*cnp+i]=ptr2[i];
      }
    } /* inner loop */

    if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
  }

  if(itno>=itmax) stop=3;

  /* restore U diagonal entries */
  for(j=mcon; j<m; ++j){
    ptr1=U + j*Usz; // set ptr1 to point to U_j
    ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
    for(i=0; i<cnp; ++i)
      ptr1[i*cnp+i]=ptr2[i];
  }

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=ea_inf;
    info[3]=dp_L2;
    for(j=mcon, tmp=DBL_MIN; j<m; ++j){
      ptr1=U + j*Usz; // set ptr1 to point to U_j
      for(i=0; i<cnp; ++i)
        if(tmp<ptr1[i*cnp+i]) tmp=ptr1[i*cnp+i];
    }
    info[4]=mu/tmp;
    info[5]=itno;
    info[6]=stop;
    info[7]=nfev;
    info[8]=njev;
    info[9]=nlss;
  }
  //sba_print_sol(n, m, p, cnp, 0, x, mnp, &idxij, rcidxs, rcsubs);
  retval=(stop!=4)?  itno : -1;
                                                               
freemem_and_return: /* NOTE: this point is also reached via a goto! */

   /* free whatever was allocated */
  free(jac); free(U);
  free(e); free(ea);  
  free(dp);
  free(rcidxs); free(rcsubs);

  free(hx); free(diagU); free(pdp);
  if(fdj_data.hxx){ // cleanup
    free(fdj_data.hxx);
    free(fdj_data.func_rcidxs);
  }

  sba_crsm_free(&idxij);

  return retval;
}


/* Bundle adjustment on structure parameters only 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 */

int sba_str_levmar_x(
    const int n,   /* number of points */
    const int m,   /* number of images */
    char *vmask,  /* visibility mask: vmask[i][j]=1 if point i visible in image j, 0 otherwise. nxm */
    double *p,    /* initial parameter vector p0: (b1, ..., bn).
                   * bi are the i-th point parameters, * size n*pnp */
    const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
    double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                   * x_ij is the projection of the i-th point on the j-th image.
                   * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                   * see vmask[i][j], max. size n*m*mnp
                   */
    const int mnp,/* number of parameters for EACH measurement; usually 2 */
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                                              /* functional relation describing measurements. Given a parameter vector p,
                                               * computes a prediction of the measurements \hat{x}. p is (n*pnp)x1,
                                               * \hat{x} is (n*m*mnp)x1, maximum
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               */
    void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                                              /* function to evaluate the sparse jacobian dX/dp.
                                               * The Jacobian is returned in jac as
                                               * (dx_11/db_1, ..., dx_1m/db_1, ..., dx_n1/db_n, ..., dx_nm/db_n), or (using HZ's notation),
                                               * jac=(B_11, ..., B_1m, ..., B_n1, ..., B_nm)
                                               * Notice that depending on idxij, some of the B_ij might be missing.
                                               * Note also that B_ij are mnp x pnp matrices and they
                                               * should be stored in jac in row-major order.
                                               * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                                               * as working memory
                                               *
                                               * If NULL, the jacobian is approximated by repetitive func calls and finite
                                               * differences. This is computationally inefficient and thus NOT recommended.
                                               */
    void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

    int itmax,         /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
    int verbose,       /* I: verbosity */
    double opts[SBA_OPTSSZ],
	                     /* I: minim. options [\mu, \epsilon1, \epsilon2]. Respectively the scale factor for initial \mu,
                        * stopping thresholds for ||J^T e||_inf, ||dp||_2 and ||e||_2
                        */
    double info[SBA_INFOSZ]
	                     /* O: information regarding the minimization. Set to NULL if don't care
                        * info[0]=||e||_2 at initial p.
                        * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                        * info[5]= # iterations,
                        * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                        *                                 2 - stopped by small dp
                        *                                 3 - stopped by itmax
                        *                                 4 - singular matrix. Restart from current p with increased mu 
                        *                                 5 - stopped by small ||e||_2
                        *                                 6 - too many attempts to increase damping. Restart with increased mu
                        * info[7]= # function evaluations
                        * info[8]= # jacobian evaluations
			                  * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                        */
)
{
register int i, j, ii, jj, k;
int nvis, nnz, retval;

/* The following are work arrays that are dynamically allocated by sba_str_levmar_x() */
double *jac;  /* work array for storing the jacobian, max. size n*m*mnp*pnp */
double *V;    /* work array for storing the V_i in the order V_1, ..., V_n, size n*pnp*pnp */

double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                 max. size n*m*mnp */
double *eb;   /* work array for storing the eb_i in the order eb_1, .. eb_n size n*pnp */

double *dp;   /* work array for storing the parameter vector updates db_1, ..., db_n, size n*pnp */

/* Of the above arrays, jac, e, are sparse and
 * V, eb, dp are dense. Sparse arrays are indexed through
 * idxij (see below), that is with the same mechanism as the input 
 * measurements vector x
 */

/* submatrices sizes */
int Bsz, Vsz,
    esz, ebsz;

register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location
                        * of B_ij in jac, etc.
                        * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                        * index k and it is used to efficiently lookup the memory locations where the non-zero
                        * blocks of a sparse matrix/vector are stored
                        */
int maxnm=(n>=m)? n:m, /* max. of (n, m) */
    *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                 column in a sparse matrix, size max(n, m) */
    *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                 sparse matrix, size max(n, m) */

/* The following variables are needed by the LM algorithm */
register int itno;  /* iteration counter */
int nsolved;
/* temporary work arrays that are dynamically allocated */
double *hx,         /* \hat{x}_i, max. size m*n*mnp */
       *diagV,      /* diagonals of V_i, size n*pnp */
       *pdp;        /* p + dp, size n*pnp */

register double mu,  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
double p_eL2, eb_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
double p_L2, dp_L2=DBL_MAX, dF, dL;
double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2], eps3_sq=opts[3]*opts[3];
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0, nlss=0;
int nobs, nvars;

struct fdj_data_x_ fdj_data;
void *jac_adata;

/* Initialization */

  /* block sizes */
  Bsz=mnp * pnp; Vsz=pnp * pnp;
  esz=mnp; ebsz=pnp;

  /* count total number of visible image points */
  for(i=nvis=0, jj=n*m; i<jj; ++i)
    nvis+=vmask[i];

  nobs=nvis*mnp;
  nvars=n*pnp;
  if(nobs<nvars){
    fprintf(stderr, "sba_str_levmar_x(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
    exit(1);
  }

  /* allocate work arrays */
  jac=(double *)emalloc(nvis*Bsz*sizeof(double));
  V=(double *)emalloc(n*Vsz*sizeof(double));
  e=(double *)emalloc(nobs*sizeof(double));
  eb=(double *)emalloc(nvars*sizeof(double));
  dp=(double *)emalloc(nvars*sizeof(double));
  rcidxs=(int *)emalloc(maxnm*sizeof(int));
  rcsubs=(int *)emalloc(maxnm*sizeof(int));


  hx=(double *)emalloc(nobs*sizeof(double));
  diagV=(double *)emalloc(nvars*sizeof(double));
  pdp=(double *)emalloc(nvars*sizeof(double));

  /* allocate & fill up the idxij structure */
  sba_crsm_alloc(&idxij, n, m, nvis);
  for(i=k=0; i<n; ++i){
    idxij.rowptr[i]=k;
    ii=i*m;
    for(j=0; j<m; ++j)
      if(vmask[ii+j]){
        idxij.val[k]=k;
        idxij.colidx[k++]=j;
      }
  }
  idxij.rowptr[n]=nvis;

  /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
  if(!fjac){
    fdj_data.func=func;
    fdj_data.cnp=0;
    fdj_data.pnp=pnp;
    fdj_data.mnp=mnp;
    fdj_data.hx=hx;
    fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
    fdj_data.func_rcidxs=(int *)emalloc(2*maxnm*sizeof(int));
    fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxnm;
    fdj_data.adata=adata;

    fjac=sba_fdjac_x;
    jac_adata=(void *)&fdj_data;
  }
  else{
    fdj_data.hxx=NULL;
    jac_adata=adata;
  }

  if(itmax==0){ /* verify jacobian */
    sba_str_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, pnp, mnp, adata, jac_adata);
    retval=0;
    goto freemem_and_return;
  }

  /* compute the error vectors e_ij in hx */
  (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
  /* compute e=x - f(p) and its L2 norm */
  for(i=0, p_eL2=0.0; i<nobs; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }

  if(verbose) printf("initial str-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
  init_p_eL2=p_eL2;

  for(itno=stop=0; itno<itmax && !stop; ++itno){
    /* Note that p, e and ||e||_2 have been updated at the previous iteration */

    /* compute derivative submatrices B_ij */
    (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

    /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
    /* V_i is symmetric, therefore its computation can be speeded up by
     * computing only the upper part and then reusing it for the lower one.
     * Recall that B_ij is mnp x pnp
     */
    /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
    /* Recall that e_ij is mnp x 1
     */
	  _dblzero(V, n*Vsz); /* clear all V_i */
	  _dblzero(eb, n*ebsz); /* clear all eb_i */
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      ptr2=eb + i*ebsz; // set ptr2 to point to eb_i

      nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
      for(j=0; j<nnz; ++j){
        /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
        ptr3=jac + idxij.val[rcidxs[j]]*Bsz;
      
        /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
        for(ii=0; ii<pnp; ++ii){
          for(jj=ii; jj<pnp; ++jj){
            for(k=0, sum=0.0; k<mnp; ++k)
              sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
            ptr1[ii*pnp+jj]+=sum;
          }

          /* copy the LOWER TRIANGULAR PART of V_i from the upper one */
          for(jj=0; jj<ii; ++jj)
            ptr1[ii*pnp+jj]=ptr1[jj*pnp+ii];
        }

        ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
        /* compute B_ij^T e_ij and add it to eb_i */
        for(ii=0; ii<pnp; ++ii){
          for(jj=0, sum=0.0; jj<mnp; ++jj)
            sum+=ptr3[jj*pnp+ii]*ptr4[jj];
          ptr2[ii]+=sum;
        }
      }
    }

    /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=0, p_L2=eb_inf=0.0; i<nvars; ++i){
      if(eb_inf < (tmp=FABS(eb[i]))) eb_inf=tmp;
      p_L2+=p[i]*p[i];
    }
    //p_L2=sqrt(p_L2);

    /* save diagonal entries so that augmentation can be later canceled.
     * Diagonal entries are in V_i
     */
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
      for(j=0; j<pnp; ++j)
        ptr2[j]=ptr1[j*pnp+j];
    }

/*
if(!(itno%100)){
  printf("Current estimate: ");
  for(i=0; i<nvars; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g\n", eb_inf, p_eL2);
}
*/

    /* check for convergence */
    if((eb_inf <= eps1)){
      dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(itno==0){
      for(i=0, tmp=DBL_MIN; i<nvars; ++i)
        if(diagV[i]>tmp) tmp=diagV[i]; /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */
    while(1){
      /* augment V */
      for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        for(j=0; j<pnp; ++j)
          ptr1[j*pnp+j]+=mu;
      }

      /* solve the linear systems V*_i db_i = eb_i to compute the db_i */
      for(i=nsolved=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
        ptr3=dp + i*pnp; // set ptr3 to point to db_i

        //nsolved+=sba_Axb_LU(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
        nsolved+=sba_Axb_Chol(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
        //nsolved+=sba_Axb_BK(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
        //nsolved+=sba_Axb_QRnoQ(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
        //nsolved+=sba_Axb_QR(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
        //nsolved+=sba_Axb_SVD(ptr1, ptr2, ptr3, pnp, 0); ++nlss;
	      //nsolved+=(sba_Axb_CG(ptr1, ptr2, ptr3, pnp, (3*pnp)/2, 1E-10, SBA_CG_JACOBI, 0) > 0); ++nlss;
      }

      if(nsolved==n){

        /* parameter vector updates are now in dp */

        /* compute p's new estimate and ||dp||^2 */
        for(i=0, dp_L2=0.0; i<nvars; ++i){
          pdp[i]=p[i] + (tmp=dp[i]);
          dp_L2+=tmp*tmp;
        }
        //dp_L2=sqrt(dp_L2);

        if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
       //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
         stop=4;
         break;
       }

        (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
        if(verbose>1)
          printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
        for(i=0, pdp_eL2=0.0; i<nobs; ++i){ /* compute ||e(pdp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pdp_eL2+=tmp*tmp;
        }

        for(i=0, dL=0.0; i<nvars; ++i)
          dL+=dp[i]*(mu*dp[i]+eb[i]);

        dF=p_eL2-pdp_eL2;

        if(verbose>1)
          printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(2.0*dF/dL-1.0);
          tmp=1.0-tmp*tmp*tmp;
          mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
          nu=2;

          for(i=0; i<nvars; ++i) /* update p's estimate */
            p[i]=pdp[i];

          for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pdp_eL2;
          break;
        }
      } /* nsolved==n */

      /* if this point is reached, either at least one linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown) */
        fprintf(stderr, "Too many failed attempts to increase the damping factor in sba_str_levmar_x()! Singular hessian matrix?\n");
        //exit(1);
        stop=6;
        break;
      }
      nu=nu2;

      /* restore V diagonal entries */
      for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
        for(j=0; j<pnp; ++j)
          ptr1[j*pnp+j]=ptr2[j];
      }
    } /* inner loop */

    if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
  }

  if(itno>=itmax) stop=3;

  /* restore V diagonal entries */
  for(i=0; i<n; ++i){
    ptr1=V + i*Vsz; // set ptr1 to point to V_i
    ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
    for(j=0; j<pnp; ++j)
      ptr1[j*pnp+j]=ptr2[j];
  }

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=eb_inf;
    info[3]=dp_L2;
    for(i=0; i<n; ++i){
      ptr1=V + i*Vsz; // set ptr1 to point to V_i
      for(j=0; j<pnp; ++j)
        if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
      }
    info[4]=mu/tmp;
    info[5]=itno;
    info[6]=stop;
    info[7]=nfev;
    info[8]=njev;
    info[9]=nlss;
  }
  //sba_print_sol(n, m, p, 0, pnp, x, mnp, &idxij, rcidxs, rcsubs);
  retval=(stop!=4)?  itno : -1;
                                                               
freemem_and_return: /* NOTE: this point is also reached via a goto! */

   /* free whatever was allocated */
  free(jac); free(V);
  free(e); free(eb);  
  free(dp);               
  free(rcidxs); free(rcsubs);

  free(hx); free(diagV); free(pdp);
  if(fdj_data.hxx){ // cleanup
    free(fdj_data.hxx);
    free(fdj_data.func_rcidxs);
  }

  sba_crsm_free(&idxij);

  return retval;
}
