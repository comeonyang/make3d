//-*-c++-*-

#ifndef _GENERALRECTIFYPLANAR_H
#define _GENERALRECTIFYPLANAR_H

#include <vector>
#include <math.h>
#include <utility>

#include "Ops.h"

namespace Rectify {

  /*! This class encapsulates the process of rectifying a pair of images using a general rectification method.
   * Order of use is:
   *    1.) Call constructor with image details and some point matches
   *    2.) If desired resample the images with resampleIms to get the actual rectified images
   *    3.) Use RectifyPoint & UnRectifyPoint to move points to and from the resampled images
   */
  class GeneralPlanarRectify {
  private:

    // These variables are for the nonlinear rectification stage
    Geometry::Homg2DPoint<double> epipole;
    // Rectification tables for nonlinear scaling in y
    typedef std::vector<std::pair<double, double> > RectifyTabT;
    typedef std::vector<double> UnRectTabT;
    RectifyTabT RectTab;
    UnRectTabT UnRectTab;

    // General Variables
    // Image sizes
    unsigned int im1xs, im1ys, im2xs, im2ys;
    // Max and min angle in the image
    double minang, maxang;
    // X region occupied by rectified images
    double minxim1, maxxim1, minxim2, maxxim2;
    // This variable indicates the images should be flipped about the x axis
    bool flipx;

    // For the compatible homographies stuff
    // The fundamental matrix
    const MultiViewGeom::FMatrix<double> &FM;
    // Homography representing the zero disparity plane
    // This transfers points from the given image to the rectify image (whichever
    // that is).
    // For the rectify image these will be the identity matrix
    MatVec::Matrix<double> Plane0Im1, Plane0InvIm1;
    MatVec::Matrix<double> Plane0Im2, Plane0InvIm2;

    // Add two angles with appropriate wrap around at 360 degrees
    double addtoangle(double angle, const double addon) const;

    // --------------------------------------------
    // Functions operating in the unrectified image
    // --------------------------------------------

    // Get the distance of the epipole from the image centre (pass image x and y sizes)
    // Returns DBL_MAX if it is an infinite epipole
    double GetEpiDist(const Geometry::Homg2DPoint<double> &e, const unsigned int xs, const unsigned int ys) const;
    // Gets an oriented epipolar line from a point (+ve scale factor only). Works for infinite epipole
    Geometry::Line2D<double> GetOrientELine(
      const Geometry::Homg2DPoint<double> &Point, const Geometry::Homg2DPoint<double> &e) const;
    // Transfer points using a homography
    Geometry::Point2D<double> TransferPointUsingH(const MatVec::Matrix<double> &H, const double &px, const double &py) const;
    // Transfer a line using a homography.
    Geometry::Line2D<double> TransferLineUsingH(const MatVec::Matrix<double> &H, const Geometry::Line2D<double> &Line) const;

    // -------------------------------------------
    // Functions operating in the nonlinear image.
    // -------------------------------------------
    // Get the angle a particular point makes with the given epipole
    double GetPointAngle(const Geometry::Point2D<double> &Point) const;
    // Get the distance of a point from an epipole.
    double GetDist(const Geometry::Point2D<double> &Point) const;
    // Get an epipolar line from an angle and the epipole
    Geometry::Line2D<double> GetEpipolarLineFromAng(const double angle) const;
    // Get the point in the nonlinear image given a distance and angle related to the epipole
    Geometry::Point2D<double> GetPointFromDistAndAng(const double dist, const double angle) const;

    // Convert a point in the nonlinear image to rectified y line (y) via look up table
    // and distance from epipole (x)
    Geometry::Point2D<double> RectPointNLin(const Geometry::Point2D<double> &Point) const;
    // UnRectify a point in the nonlinear image. Converts from y line via look up table
    // and distance from epipole (x) to a point in the nonlinear image
    Geometry::Point2D<double> UnRectPointNLin(const Geometry::Point2D<double> &Point) const;

    // ------------------------------------------------
    // Major functions that actualy control the process
    // ------------------------------------------------

    // Sets up the rectify table bounds
    void SetUpTablesAndXBounds(void);

    // Gets the maximal epipolar lines for the image with the given epipole and image size. Transfers these to nonlinear image using the given homography H. H need not be oriented. Uses the epipole, epipole which must be valid before calling.
    void GetImAngleRange(const unsigned int xs, const unsigned int ys, double &minang, double &maxang, const MatVec::Matrix<double> &H) const;

    //! mask of bits to & with your data to get the red channel
    static const unsigned int redmask=255;
    //! mask of bits to & with your data to get the green channel (don't forget to shift right 8 bits afterwards)
    static const unsigned int greenmask=255<<8;
    //! mask of bits to & with your data to get the blue channel (don't forget to shift right 16 bits afterwards)
    static const unsigned int bluemask=255<<16;

    // Get a sub-pixel co-ordinate in an image, using bilinear interpolation to determine a value
    template <class ImageT, class FracT>
      typename ImageT::value_type subpix_Bilinear(const ImageT &image, FracT x, FracT y) const {

        typedef typename ImageT::value_type T;

        int ix = int(x);		// integral parts
        int iy = int(y);
        FracT fx = x - ix;			// fractional parts
        FracT fy = y - iy;
        T c00 = image(ix, iy); // 4 corners of rectangle
        T c01 = image(ix, iy+1);
        T c10 = image(ix+1, iy);
        T c11 = image(ix+1, iy+1);

        if (!image.isColour()) {
          // interpolate along y
          FracT c0 = c00 + fy * (c01 - c00);
          FracT c1 = c10 + fy * (c11 - c10);
          // interpolate along x
          return (T)(c0  + fx * (c1 - c0));
        }
        else
        {
          // Do red, green and blue seperately
          FracT c00_, c01_, c11_, c10_, c0, c1;
          T rval;

          // Red
          c00_=(c00 & redmask); c01_=(c01 & redmask);
          c11_=(c11 & redmask); c10_=(c10 & redmask);
          c0 = c00_ + fy * (c01_ - c00_); 
          c1 = c10_ + fy * (c11_ - c10_);
          rval=(T)(c0  + fx * (c1 - c0));

          // Green
          c00_=(c00 & greenmask)>>8; c01_=(c01 & greenmask)>>8;
          c11_=(c11 & greenmask)>>8; c10_=(c10 & greenmask)>>8;
          c0 = c00_ + fy * (c01_ - c00_); 
          c1 = c10_ + fy * (c11_ - c10_);
          rval=rval | ((T)(c0  + fx * (c1 - c0)))<<8;

          // Blue
          c00_=(c00 & bluemask)>>16; c01_=(c01 & bluemask)>>16;
          c11_=(c11 & bluemask)>>16; c10_=(c10 & bluemask)>>16;
          c0 = c00_ + fy * (c01_ - c00_); 
          c1 = c10_ + fy * (c11_ - c10_);
          rval=rval | ((T)(c0  + fx * (c1 - c0)))<<16;

          return rval;
        }
      }

      /// Rectifies the image in ImIn to ImOut. minx, maxx give dimensions of rectified image. Uses rectify table to give height of the image.
      template <class ImageT1, class ImageT2>
        void ResampleAndRectify(
        const ImageT1 &InIm, ImageT2 &OutIm,
        const double minx, const double maxx, const bool flip,
        const MatVec::Matrix<double> &Trans,
        const typename ImageT2::value_type bound) {

          // Original unrectified image dimensions
          size_t inxs=InIm.xsize(), inys=InIm.ysize();
          // Output image dimensions
          size_t  outxs=(unsigned int)fabs(maxx-minx), outys=(unsigned int)UnRectTab.size();

          OutIm.resize(outxs, outys, InIm.isColour());

          Geometry::Point2D<double> start, end;
          double dx,dy,norm,inx,iny;
          unsigned int outx;
          Geometry::Point2D<double> TPoint;
          for (unsigned int outy=0; outy<outys; ++outy)
          {
            // Get start and end of this epipolar line in images
            start=GetPointFromDistAndAng(minx, UnRectTab[outy]);
            end=GetPointFromDistAndAng(maxx, UnRectTab[outy]);

            // If flipping the image swap start and end
            if (flip) {
              Geometry::Point2D<double> temp(start);
              start=end;
              end=temp;
            }

            // Now, work from start to end according to the gradient.
            // Normalise the gradient to have eucl. norm 1 - hence go 1 pixel at a time.
            dx=end.x()-start.x(); dy=end.y()-start.y();
            norm=1.0/hypot(dx,dy);
            dx=dx*norm; dy=dy*norm;
            inx=start.x(); iny=start.y();
            for (outx=0;outx<outxs;++outx, inx+=dx, iny+=dy)
            {
              TPoint.x()=inx; TPoint.y()=iny;
              TPoint=TransferPointUsingH(Trans, TPoint.x(), TPoint.y());
              // -2 because of interpolation
              const double tx=TPoint.x(), ty=TPoint.y();
              if (tx>=0 && tx<=inxs-2 && ty>=0 && ty<=inys-2)
                OutIm(outx, outy)=subpix_Bilinear(InIm, tx, ty);
              else
                OutIm(outx, outy)=bound;
            }
          }
        }

  public:

    //-------------------------------------
    // Public functions for use by the user
    //-------------------------------------
    /** Rectify any pair of images sharing common features using a general planar technique. To actually rectify a pair of images you'll want to call resampleIms afterwards. There is a seperation so you can recreate this object when you just want to rectify or unrectify a few points.

    {\bf On Entry:}
    \begin{itemize}
    \item{\em im1xs, im1ys, im2xs, im2ys: } The dimensions of the two images to rectify
    \item{\em FM: } The fundamental matrix for images 1,2.
    \item{\em Matches: } The point matches used to calcuate the fundamental matrix. In fact any set of point matches that are nicely distributed, outlier free and fit the epipolar geometry will do - but who gives a shit. These will be hartley sturm corrected so don't feel the need to do that
    \item{\em NLinY: } A bool, which if set to true indicates to use nonlinear scaling on the y axis of the rectified images. This nonlinear scaling is such that there will be no pixel compression and the rectified image will be as small as possible. If there are infinite epipoles this method will not be used regardless (isn't necessary for infinite epipoles anyway).
    \end{itemize}

    {\bf On Exit:}
    \begin{itemize}
    \item{\em OutIm1, OutIm2: } 2 Output images. Can not be the same as Im1, Im2. There is no need to set them up to any size or bits per pixel.
    \end{itemize}

    {\bf Errors:}
    \begin{itemize}
    \item{\em 1:} Occurs if you try to rectify two non-overlapping images. Although they can theoreticaly be rectified - there is no overlapping section, so what is the point.
    \end{itemize}
    */
    GeneralPlanarRectify(const unsigned int im1xs, const unsigned int im1ys,
      const unsigned int im2xs, const unsigned int im2ys,
      const MultiViewGeom::FMatrix<double> &FM,
      const std::vector<std::pair<Geometry::Point2D<double>, Geometry::Point2D<double> > > &Matches);

    ~GeneralPlanarRectify();

    // Call to resample the images
    // Supply input and output images as well as the value used to indicate there is no valid image (i.e. for the boundary
    // outside non-square images). I like to use int images and -1 - but use whatever you like.
    template <class ImageT, class ImageT2>
      void resampleIms(const ImageT &Im1, const ImageT &Im2,
      ImageT2 &OutIm1, ImageT2 &OutIm2,
      const typename ImageT2::value_type bound) {
        // Output images
        ResampleAndRectify(Im1, OutIm1, minxim1, maxxim1, flipx, Plane0InvIm1, bound);
        ResampleAndRectify(Im2, OutIm2, minxim2, maxxim2, flipx, Plane0InvIm2, bound);
      }

      // Rectify the supplied point (x,y). Returns the point in x,y
      void RectifyPointIm1(double &x, double &y) const;
      void RectifyPointIm2(double &x, double &y) const;
      void UnRectifyPointIm1(double &x, double &y) const;
      void UnRectifyPointIm2(double &x, double &y) const;

      // Can be used to unrectify a set of points all on the same epipolar line.
      // The epipolar line is given by the y co-ordinate of the rectified line, and the
      // input iterator dereferences to give numbers representing x co-ordinates.
      // RectPoints is added to (not erased).
      // Sorry no proper use of iterators - couldn't stand cluttering up the .h file
      void UnRectifyPointsIm1(std::vector<double> x, const double y, std::vector<Geometry::Point2D<double> > &RectPoints) const;
      void UnRectifyPointsIm2(std::vector<double> x, const double y, std::vector<Geometry::Point2D<double> > &RectPoints) const;

  };

};


#endif
