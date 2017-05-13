//-*-c++-*-

#ifndef _RECTIFY_EXTRAS_H
#define _RECTIFY_EXTRAS_H

#include <assert.h>

// This files contains all the code I've had to cannabilise out of my home made libraries.
// Not a neat job - really most of these datatypes are fundamentally unnecessary and should be replaced by simple 1D arrays
// But I can't be bothered.

/*******************************/
/* Images and matrices */



/*

// EXAMPLE!
// If you want to use your image type directly then you can use an adaptor as below.
// This example shows adaption of a user image type MyImage that has equivalent functions
// recreate, isrgb, getpixel, xdim, ydim

namespace Image {

  class ImageAdaptor {
  protected:

    MyImage &image;

  public:

    typedef unsigned int value_type;

    ImageAdaptor(MyImage &im) : image(im) {}

    void resize(size_t w, size_t h, bool iscol) {
      image.recreate(w,h);
    }

    inline unsigned int& operator()(size_t x, size_t y) { return image.getpix(x,y); }
    inline const unsigned int& operator()(size_t x, size_t y) const { return image.getpix(x,y); }

    bool isColour() const { return image.isrgb(); }

    inline size_t xsize() const { return image.xdim(); }
    inline size_t ysize() const { return image.ydim(); }
  }
}
*/




namespace Image {

/*
  This class demonstrates all the functionality an image class must offer if it is to be supplied to the rectification routines
  It works as a standalone very very basic image class too.
  If you have an image class then simply use this as a wrapper (i.e. make the constructor take your image class and redirect
  all other functions to your class). See above for an example
 */
template <class T=unsigned int>
class Image {
protected:

  T *data;
  size_t xs, ys;
  bool iscolour;

public:


  /// The type of object stored in the container
  typedef T value_type;

  Image() : data(NULL), xs(0), ys(0), iscolour(false) { }
  Image(size_t w, size_t h, bool colour) : data(new T[w*h]), xs(w), ys(h), iscolour(colour) {
    assert(w>0 && h>0);
  }

  ~Image() {
    if (data!=NULL)
      delete[] data;
  }
  void resize(size_t w, size_t h, bool iscol) {
    assert(w>0 && h>0);
    iscolour=iscol;
    if (w!=xs || h!=ys) {
      if (data!=NULL)
	delete[] data;

      xs=w; ys=h;
      data=new T[xs*ys];
    }
  }

  inline T& operator()(size_t x, size_t y) {
    return data[x+y*xs];
  }
  inline const T& operator()(size_t x, size_t y) const {
    return data[x+y*xs];
  }

  bool isColour() const { return iscolour; }

  inline size_t xsize() const { return xs; }
  inline size_t ysize() const { return ys; }

};

}

namespace MatVec {

  template <class T>
  class Matrix {
  protected:

    T *data;
    size_t nrows, ncols, length;
    bool iscolour;

    void copy(const Matrix<T> &Other) {
      resize(Other.nrows, Other.ncols, Other.iscolour);
      if (length!=0) {
        T *begin=Other.data, *end=Other.data+length, *out=data;
        for (;begin<end;)
          *out++=*begin++;
      }
    }

  public:

    /// The type of object stored in the container
    typedef T value_type;

    Matrix() : data(NULL), nrows(0), ncols(0), length(0), iscolour(false) { }
    Matrix(size_t nrows_, size_t ncols_) : data(new T[nrows_*ncols_]), nrows(nrows_), ncols(ncols_), length(nrows_*ncols) {
      assert(nrows>0 && ncols>0);
    }
    /// Copy constructor. Performs a deep copy
    Matrix(const Matrix<T> &Other) : data(NULL), nrows(0), ncols(0), length(0) { copy(Other); }
    /// Allow = assignment
    Matrix<T>& operator=(const Matrix<T> &Other) { copy(Other); return *this;}

    ~Matrix() {
      if (data!=NULL)
        delete[] data;
    }

    void resize(size_t w, size_t h, bool colour) {
      iscolour=colour;
      assert(w>0 && h>0);
      if (w!=nrows || h!=ncols) {
        if (data!=NULL)
          delete[] data;

        nrows=w; ncols=h; length=nrows*ncols;
        data=new T[nrows*ncols];
      }
    }

    inline T& operator()(size_t x, size_t y) {
      assert(x<nrows && y<ncols);
      return data[x+y*nrows];
    }
    inline const T& operator()(size_t x, size_t y) const {
      assert(x<nrows && y<ncols);
      return data[x+y*nrows];
    }

    inline size_t num_rows() const { return nrows; }
    inline size_t num_cols() const { return ncols; }
  };

}

namespace MultiViewGeom {

  template <class T>
  class FMatrix: public MatVec::Matrix<T> {
  public:

    FMatrix() : MatVec::Matrix<T>(3,3) {}
  };

}




























/*******************************/
/* Geometry */


namespace Geometry {

/*! \brief Very very simple point classes
*/

template <class T>
class Point2D {
private:
  T vals[2];

public:
  typedef T value_type;

  Point2D() {}
  Point2D(T x, T y) { vals[0]=x; vals[1]=y; }
  Point2D(const Point2D<T> &P) { vals[0]=P.vals[0]; vals[1]=P.vals[1]; }
  /// Specialised = assignment (for efficiency and to stop cfront errors).
  Point2D &operator=(const Point2D<T> &P) {
    vals[0]=P.vals[0]; vals[1]=P.vals[1];
    return *this;
  }

  /// Get direct access to x
  inline T& x() { return vals[0]; }
  /// Get direct access to x
  inline T& y() { return vals[1]; }
  /// Can't overwrite access to x
  inline const T& x() const { return vals[0]; }
  /// Can't overwrite access to x
  inline const T& y() const { return vals[1]; }

  inline size_t size() const { return 2; }
  inline T& operator[](size_t n) { assert(n<2); return vals[n]; }
  inline const T& operator[](size_t n) const { assert(n<2); return vals[n]; }
};

template <class T>
class Homg2DPoint {
private:
  T vals[3];

public:
  typedef T value_type;

  Homg2DPoint() { vals[0]=0; vals[1]=0; vals[2]=0;}
  Homg2DPoint(T x, T y, T w) { vals[0]=x; vals[1]=y; vals[2]=w;}
  Homg2DPoint(const Homg2DPoint<T> &P) { vals[0]=P.vals[0]; vals[1]=P.vals[1]; vals[2]=P.vals[2]; }
  /// Specialised = assignment (for efficiency and to stop cfront errors).
  Homg2DPoint &operator=(const Homg2DPoint<T> &P) { vals[0]=P.vals[0]; vals[1]=P.vals[1]; vals[2]=P.vals[2]; return *this; }

  inline T& x() { return vals[0]; }
  inline T& y() { return vals[1]; }
  inline T& w() { return vals[2]; }
  inline const T& x() const { return vals[0]; }
  inline const T& y() const { return vals[1]; }
  inline const T& w() const { return vals[2]; }

  inline T& operator[](size_t n) { assert(n<3); return vals[n]; }
  inline const T& operator[](size_t n) const { assert(n<3); return vals[n]; }

  inline size_t size() const { return 3; }
};

template <class T>
class Line2D {
private:
  T vals[3];

public:
  typedef T value_type;

  Line2D() { }
  Line2D(T a, T b, T c) { vals[0]=a; vals[1]=b; vals[2]=c;}
  Line2D(const Line2D<T> &P) { vals[0]=P.vals[0]; vals[1]=P.vals[1]; vals[2]=P.vals[2]; }
  /// Specialised = assignment (for efficiency and to stop cfront errors).
  Line2D &operator=(const Line2D<T> &P) { vals[0]=P.vals[0]; vals[1]=P.vals[1]; vals[2]=P.vals[2]; return *this; }

  inline T& a() { return vals[0]; }
  inline T& b() { return vals[1]; }
  inline T& c() { return vals[2]; }
  inline const T& a() const { return vals[0]; }
  inline const T& b() const { return vals[1]; }
  inline const T& c() const { return vals[2]; }

  inline size_t size() const { return 3; }
  inline T& operator[](size_t n) { assert(n<3); return vals[n]; }
  inline const T& operator[](size_t n) const { assert(n<3); return vals[n]; }
};

}


























/*******************************/
/*! Calculate homographies that are consistant with a particular fundamental matrix. These use extra matches to select a homography from the 3 parameter family that produce a mapping that is consistant with the fundamental matrix (i.e. homographies that map points on epipolar lines to points on epipolar lines).
 */
namespace HomogConsistF {

  void RANSAC(const std::vector<std::pair<Geometry::Point2D<double>, Geometry::Point2D<double> > > &Matches,
	      std::vector<std::pair<Geometry::Point2D<double>, Geometry::Point2D<double> > > &OutlierFree,
	      MatVec::Matrix<double> &H, const MultiViewGeom::FMatrix<double> &FM);

}


// ******************
// NUMERICAL ROUTINES
// ******************

namespace MatrixOps {

  extern void EigenValVec_Symmetric(const MatVec::Matrix<double> &A, double *EigenVals, MatVec::Matrix<double> &EigenVecs);
  extern bool inverse(const MatVec::Matrix<double> &A, MatVec::Matrix<double> &AInv);
}

namespace LAPACK {

  // allow version with vectors. NOTE: B & X Can safely be the same vector
  // Overwrites A
  extern bool LeastSquaresRankDeff(MatVec::Matrix<double>& A, const double *B, double *X);

}

/** This function can be called to get the median from the given iterator range. Assumes is unsorted. */
namespace MiscMath {

  extern double VectorMedian(double *VecBegin, double *VecEnd);

}

#endif
