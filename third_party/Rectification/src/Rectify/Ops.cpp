#include <math.h>
#include <string>
#include <vector>
#include <float.h>
#include <algorithm>
#ifndef _WIN32
#include <ext/algorithm>
#endif
#include <iostream>

#include "Ops.h"

//#define DEBUG_RANSAC

using namespace MatVec;
using namespace std;
using namespace Geometry;
using namespace MultiViewGeom;
typedef vector<pair<Point2D<double>, Point2D<double> > > MatchVec2Im;

namespace HomogConsistF {

  // Linear version for use in RANSAC
  template <class T, class MatchT>
    bool Linear(Matrix<T> &H, MatchT &Matches, const MultiViewGeom::FMatrix<T> &FM) {
      Matrix<double> EigenVecs(3,3), FT(3,3);
      // Do FM^{T}*FM
      double sum;
      // TransposeT()
      for (unsigned int i=0;i<3;i++)
        for (unsigned int k=0;k<3;k++)
        {
          sum=FM(i,0)*FM(k,0);
          for (unsigned int j=1;j<3;j++)
            sum+=FM(i,j)*FM(k,j);

          FT(i,k)=sum;
        }
        double EigenVals[3];
        MatrixOps::EigenValVec_Symmetric(FT, EigenVals, EigenVecs);
        Homg2DPoint<double> e2(EigenVecs(0,0), EigenVecs(1,0), EigenVecs(2,0));
        double norm=1.0/sqrt(e2.x()+e2.y()+e2.w());
        e2.x()*=norm; e2.y()*=norm; e2.w()*=norm;

        const T e2x=e2.x(), e2y=e2.y(), e2w=e2.w();
        // Helps numerical stability in some near degenerate cases to at least have a value in a1,a2,a3
        T a1=0.3, a2=0.3, a3=0.3;
        H(0,0)=e2y*FM(2,0)-e2w*FM(1,0)+e2x*a1;
        H(0,1)=e2y*FM(2,1)-e2w*FM(1,1)+e2x*a2;
        H(0,2)=e2y*FM(2,2)-e2w*FM(1,2)+e2x*a3;

        H(1,0)=e2w*FM(0,0)-e2x*FM(2,0)+e2y*a1;
        H(1,1)=e2w*FM(0,1)-e2x*FM(2,1)+e2y*a2;
        H(1,2)=e2w*FM(0,2)-e2x*FM(2,2)+e2y*a3;

        H(2,0)=e2x*FM(1,0)-e2y*FM(0,0)+e2w*a1;
        H(2,1)=e2x*FM(1,1)-e2y*FM(0,1)+e2w*a2;
        H(2,2)=e2x*FM(1,2)-e2y*FM(0,2)+e2w*a3;

        // Build a constraint matrix
        Matrix<T> A(Matches.size(), 3);
        double *b=new double[Matches.size()];

        for (unsigned int i=0; i<Matches.size(); ++i) {
          const T ux=Matches[i].first.x(), uy=Matches[i].first.y();
          const T u2x=Matches[i].second.x(), u2y=Matches[i].second.y();
          const T u2w=1.0;

          const T theta2=-atan((u2w*e2y-e2w*u2y)/(u2w*e2x-e2w*u2x));

          // Constraints on x
          A(i,0)=(-cos(theta2)*u2w*e2x+sin(theta2)*u2w*e2y-(-cos(theta2)*u2x+sin(theta2)*u2y)*e2w)*ux;
          A(i,1)=(-cos(theta2)*u2w*e2x+sin(theta2)*u2w*e2y-(-cos(theta2)*u2x+sin(theta2)*u2y)*e2w)*uy;
          A(i,2)=(-cos(theta2)*u2w*e2x+sin(theta2)*u2w*e2y-(-cos(theta2)*u2x+sin(theta2)*u2y)*e2w);
          b[i]=-((cos(theta2)*u2w*H(0,0)-sin(theta2)*u2w*H(1,0)+(-cos(theta2)*u2x+sin(theta2)*u2y)*H(2,0))*ux+(-cos(theta2)*u2x+sin(theta2)*u2y)*H(2,2)+(cos(theta2)*u2w*H(0,1)-sin(theta2)*u2w*H(1,1)+(-cos(theta2)*u2x+sin(theta2)*u2y)*H(2,1))*uy+cos(theta2)*u2w*H(0,2)-sin(theta2)*u2w*H(1,2));

        }

        double x[3];
        bool ok=LAPACK::LeastSquaresRankDeff(A,b,x);

        delete[] b;

        a1=x[0]; a2=x[1]; a3=x[2];

        H(0,0)=H(0,0)-e2x*a1;
        H(0,1)=H(0,1)-e2x*a2;
        H(0,2)=H(0,2)-e2x*a3;

        H(1,0)=H(1,0)-e2y*a1;
        H(1,1)=H(1,1)-e2y*a2;
        H(1,2)=H(1,2)-e2y*a3;

        H(2,0)=H(2,0)-e2w*a1;
        H(2,1)=H(2,1)-e2w*a2;
        H(2,2)=H(2,2)-e2w*a3;

        return ok;
    }












    double GetResidual_RANSAC(const pair<Point2D<double>, Point2D<double> > &Match,
      const Matrix<double> &H) {

        const double px=Match.first[0], py=Match.first[1], pw=1.0;
        double tx=H(0,0)*px+H(0,1)*py+H(0,2)*pw;
        double ty=H(1,0)*px+H(1,1)*py+H(1,2)*pw;
        double tz=H(2,0)*px+H(2,1)*py+H(2,2)*pw;
        Point2D<double> NP(tx/tz,ty/tz);
        double xdiff=NP.x()-Match.second.x(), ydiff=NP.y()-Match.second.y();
        return xdiff*xdiff+ydiff*ydiff;
      }

      // The function itself
      void RANSAC(const vector<pair<Point2D<double>, Point2D<double> > > &Matches,
        vector<pair<Point2D<double>, Point2D<double> > > &OutlierFree,
        Matrix<double> &ReturnH, const FMatrix<double> &FM) {

          // Copy out of the iterator and into a vector so can normalise without destroying callers copy.
          unsigned int i;
          vector<Matrix<double> > Hs;
          // I know - I know - this is just in this simplified version
          OutlierFree.clear();

          vector<pair<Point2D<double>, Point2D<double> > >::size_type numm=Matches.size();

          if (numm<10) {
            cerr << "Planar homography RANSAC estimate (consistant with fundamental matrix) must be given at least 10 point matches to work!. So, go away and don't come back until you've got some more!\n" << endl;
            exit(1);
          }

          unsigned int minsample=10;
          // LMedS style calculation of number of sub-samples
          unsigned int NumberSubSamples=400;

#ifdef DEBUG_RANSAC
          cerr << "Calculating Planar Homography matrix using LMedS:" << endl;
          cerr << "Number matches=" << numm << endl;
          cerr << "Using " << NumberSubSamples << " 4 point subsamples" << endl;
#endif

          // Pre-allocate match vector for the four point algorithm
          vector<pair<Point2D<double>, Point2D<double> > > PassThrough(minsample); 

          // So, run the 3 point algorithm on NumberSubSamples subsamples
          for (i=0;i<NumberSubSamples;i++)
          {
            // They always insist on removing the most useful functions in the windows STL. B******S
#ifndef _WIN32
            // Randomly select minsample points
            random_sample_n(Matches.begin(), Matches.end(), PassThrough.begin(), minsample);
#else
            // This version is a bit naff since duplicates can occur
            unsigned int picked=0;
            while (picked<minsample) {
              int off=rand()%(int)numm;
              PassThrough[picked++]=Matches[off];
            }
#endif

            // Run the minimal algorithm
            Matrix<double> H(3,3);
            if (Linear(H, PassThrough, FM))
              Hs.push_back(H);
          }

          // Now, need to determine the median of the squared residuals, for each image
          double *residuals=new double[numm], *bestr=new double[numm];
          double minhuber=DBL_MAX, bestmaxerr=0.0;
          double med;
          int residual;
          vector<Matrix<double> >::const_iterator Hs_iter;
          MatchVec2Im::const_iterator iter;
          Matrix<double> BestH(3,3);

          for (Hs_iter=Hs.begin(); Hs_iter!=Hs.end(); ++Hs_iter)
          {
            const Matrix<double> &H=*Hs_iter;
            residual=0;

            // First run through all matches and evaluate the residual
            for (iter=Matches.begin();iter!=Matches.end();++iter)
              // evaluate d(m, Hm')^2 + d(m', H^T m)^2
              residuals[residual++]=GetResidual_RANSAC(*iter, H);

            // Now, find the median of the results.
            if (residual>0)
            {
              med=MiscMath::VectorMedian(residuals, residuals+residual);

              // Get confidence limit for Huber function
              double robuststddev=1.4826*(1.0+(5.0/(numm-minsample)))*sqrt(med);
              double maxerr=3.84*robuststddev*robuststddev;

              // Now, go through residuals and sum
              double hubersum=0.0;
              for (int i=0;i<residual;++i) {
                if (residuals[i]<maxerr)
                  hubersum+=residuals[i];
                else
                  hubersum+=maxerr;
              }

              if (hubersum<minhuber)
              {
                BestH=H; minhuber=hubersum; bestmaxerr=maxerr;
                copy(residuals, residuals+residual, bestr);
              }
            }
          }
          Hs.clear();

          // Now, discard all false matches
#ifdef DEBUG_RANSAC
          cerr << "Outlier threshold (95%) maxerr=" << bestmaxerr << endl;
#endif

          // Work through all matches and discard if residual squared is less than just
          // calculated robust value.
          unsigned int off, outliers=0;
          for (off=0;off<numm;++off)
          {
            if (bestr[off]>bestmaxerr)
              ++outliers;
            else
              OutlierFree.push_back(Matches[off]);
          }

          ReturnH=BestH;

#ifdef DEBUG_RANSAC
          cerr << endl << "Completed - found " << outliers << " outliers in " << numm << " points" << endl;
#endif

          delete[] residuals;
          delete[] bestr;
        }

#undef DEBUG_RANSAC

}

namespace MatrixOps {

  extern "C" {
    void dspevd_(char *, char *, int *, double *, double *, double *, int *, double *, int *, int *, int *, int *);
  }

  /**
  Finds the eigenvalues and eigenvectors of a symmetric matrix using LAPACK routines. Returns eigenvalues and eigenvectors in ascending order (lowest first).
  {\bf Errors:}
  \begin{itemize}
  \item {\em 1:} Error in EigenValVec_Symmetric: Supplied matrix is not symmetric - This is returned whenever the supplied matrix is not symmetric.
  \item {\em 2:} Error in EigenValVec_Symmetric: Eigenvalue and Eigenvector arrays are not the correct size - If the returning eigenvalues/eigenvector arrays are incorrectly sized this is returned..
  */
  void EigenValVec_Symmetric(const Matrix<double> &A, double *EigenVals, Matrix<double> &EigenVecs)
  {
    // First, check if it is symmetric
    unsigned int i,j;
    size_t N=A.num_rows();
    bool Symm=true;

    if (N!=A.num_cols())
      Symm=false;
    else
    {
      for (i=0;i<N-1;i++)
        for (j=0;j<i;j++)
          if (A(i,j)!=A(j,i))
            Symm=false;
    }

    if (Symm==false) {
      throw 1;
    }

    double *A_Packed=new double[(N*(N+1))/2];
    unsigned int indx=0;
    // Put into packed storage
    for (i=0;i<N;i++)
      for (j=0;j<=i;j++)
        A_Packed[indx++]=A(i,j);

    int rinfo=0;
    int LWork=(int)(1+6*N+N*N);
    double *Work=new double[LWork];
    int LIWork=(int)(3+5*N);
    int *IWork=new int[LIWork];
    int N_=(int)N; // Has to be an int for passing to fortran code

    EigenVecs.resize(N,N,false);

    dspevd_("V", "U", &N_, A_Packed, &EigenVals[0], &EigenVecs(0,0), &N_, Work, &LWork, IWork, &LIWork, &rinfo);

    delete[] A_Packed;
    delete[] Work;
    delete[] IWork;

    assert(rinfo==0);
    return;

  }

  // The LAPACK library functions
  extern "C" {
    int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
    int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  }

  bool inverse(const MatVec::Matrix<double> &A, MatVec::Matrix<double> &AInv) {

    int N = (int)A.num_rows();
    if (A.num_rows()!=A.num_cols() || N<1)
      return false;

    // Going to use AInv as our workspace
    AInv=A;

    // Set up storage
    int *ipiv = new int[N];
    int info;

    // Call for the LU factorisation
    dgetrf_(&N, &N, &(AInv(0,0)), &N, ipiv, &info);

    if (info>0) {
      std::cerr << "info=" << info << std::endl;
      delete[] ipiv;
      return false;
    }

    // Set up workspace
    int lwork = N;
    double *workd = new double[lwork];

    // Now, call for inverse
    dgetri_(&N, &(AInv(0,0)), &N, ipiv, workd, &lwork, &info);

    delete[] ipiv;
    delete[] workd;

    // If got a singular matrix
    if (info>0)
      return false;

    return true;
  }


}




namespace MiscMath {

  double selectnth_heap(unsigned long m, double *ArrBegin, double *ArrEnd) {
    typedef double DataT;

    unsigned long i,j,k,n=(unsigned long)(distance(ArrBegin, ArrEnd));
    DataT swap;
    m=n-m+1; // Since want to find nth lowest value, not nth highest value

    //  const DataT *arr=ArrBegin-1;
    DataT *heap=new DataT[m+2];

    copy(ArrBegin, ArrBegin+m, heap);
    sort(heap, heap+m);
    --heap;

    for (i=m+1;i<=n;++i) {
      if (ArrBegin[i-1]>heap[1]) {
        heap[1]=ArrBegin[i-1];
        for (j=1;;) {
          k=j << 1;
          if (k>m) break;
          if (k!=m && heap[k] > heap[k+1]) ++k;
          if (heap[j] <= heap[k]) break;
          swap=heap[k];
          heap[k]=heap[j];
          heap[j]=swap;
          j=k;
        }
      }
    }

    DataT rval=heap[1];
    ++heap;
    delete[] heap;

    return rval;
  }

  double VectorMedian(double *VecBegin, double *VecEnd) {
    if (VecBegin==VecEnd)
      return 0.0;

    double median;

    unsigned int size=(unsigned int)(distance(VecBegin, VecEnd));

    if (size % 2 == 0)
    {
      // Don't use pedants if enough items in the iterator to get away with it
      if (size>100)
        median=MiscMath::selectnth_heap(size/2, VecBegin, VecEnd);
      else
        median=
        (MiscMath::selectnth_heap(size/2, VecBegin, VecEnd)+
        MiscMath::selectnth_heap((size/2)+1, VecBegin, VecEnd))/2.0;
    }
    else
      // Odd case
      median=MiscMath::selectnth_heap(((size+1)/2)-1, VecBegin, VecEnd);

    return median;
  }

}

namespace LAPACK {

  extern "C" {
    int dgelsy_(int *, int *, int *, double *, int *, double *, int *, int *, double *, int *, double *, int *, int *);
  }

  bool LeastSquaresRankDeff(Matrix<double>& A, const double *B, double *X) {
    int m=(int)A.num_rows(), n=(int)A.num_cols();
    int info=1, one=1;
    int *jpvt=new int[n];
    int rank;
    double rcond=1e-14;

    int mn=m;
    if (n<m)
      mn=n;
    int lwork=mn+3*n+1;
    if (2*mn+1>lwork)
      lwork=2*mn+1;
    double *Work=new double[lwork];

    // copy B into temp holder since X is wrong dimensions
    double *Temp=new double[m];
    for (int i=0;i<m;++i)
      Temp[i]=B[i];

    dgelsy_(&m, &n, &one, &(A(0,0)), &m, &Temp[0], &m, &jpvt[0], &rcond, &rank, &Work[0], &lwork, &info);

    if (info<0) {
      //      cerr << "Error from linear least square solver (LAPACK routines)" << endl;
      return false;
    }

    for (int i=0;i<n;++i)
      X[i]=Temp[i];

    delete[] jpvt;
    delete[] Work;
    delete[] Temp;
    return true;
  }

}
