// needed under WIN32 to get M_PI defined because M_PI is clearly such an unusual constant.
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <float.h>
#include <iostream>

#include "General2Im.h"

//#define DEBUG_GENERALRECTIFY

using namespace MatVec;
using namespace MultiViewGeom;
using namespace Geometry;
using namespace std;

namespace Rectify {

  void GetEpipoles(const FMatrix<double> &FM, Homg2DPoint<double> &e1, Homg2DPoint<double> &e2) {
    double EigenVals[3];
    EigenVals[0]=0.0; // Just to stop compiler warnings - it IS used.
    Matrix<double> EigenVecs(3,3), FT(3,3);

    // Do FM*FM^{T}
    // TTranspose
    double sum;
    for (unsigned int i=0;i<3;++i) {
      for (unsigned int k=0;k<=i;++k)
      {
        sum=FM(0,i)*FM(0,k);
        for (unsigned int j=1;j<3;++j)
          sum+=FM(j,i)*FM(j,k);
        FT(i,k)=sum;
        FT(k,i)=sum;
      }
    }

    MatrixOps::EigenValVec_Symmetric(FT, EigenVals, EigenVecs);
    e1[0]=EigenVecs(0,0);
    e1[1]=EigenVecs(1,0);
    e1[2]=EigenVecs(2,0);

    // Do FM^{T}*FM
    // TransposeT()
    for (unsigned int i=0;i<3;i++) {
      for (unsigned int k=0;k<3;k++)
      {
        sum=FM(i,0)*FM(k,0);
        for (unsigned int j=1;j<3;j++)
          sum+=FM(i,j)*FM(k,j);

        FT(i,k)=sum;
      }
    }


    MatrixOps::EigenValVec_Symmetric(FT, EigenVals, EigenVecs);
    e2[0]=EigenVecs(0,0);
    e2[1]=EigenVecs(1,0);
    e2[2]=EigenVecs(2,0);
  }


  // -----------------------------------------
  // Functions extracted from my geometry code
  // to make this file less dependent
  // These are a bit naff!
  // -----------------------------------------
  template <class T>
  inline T MIN(T a, T b) {
    if (a<b)
      return a;
    return b;
  }

  template <class T>
  inline T MAX(T a, T b) {
    if (a>b)
      return a;
    return b;
  }

  template <class T>
  inline short SGN(T a) {
    if (a<0.0)
      return -1;
    else
      return 1;
  }

  template <class VectorT1, class VectorT2, class VectorT3>
  inline void CrossProduct(const VectorT1 &a, const VectorT2 &b, VectorT3 &result) {
    assert(b.size()==3 && a.size()==3 && result.size()==3);
    result[0]=a[1]*b[2]-a[2]*b[1];
    result[1]=a[2]*b[0]-a[0]*b[2];
    result[2]=a[0]*b[1]-a[1]*b[0];
  }

  /// Get the distance between 2 2D points.
  template <class Point2DT, class Point2DT2>
  static inline typename Point2DT::value_type distance_points(const Point2DT &P1, const Point2DT2 &P2) {
    return hypot(P1[0]-P2[0], P1[1]-P2[1]);
  }

  /// Test if parallel with another line.
  template <class Line2DT1, class Line2DT2, class T>
  static bool Parallel(const Line2DT1 &Line1, const Line2DT2 &Line2, const T tol) {
    T d1x, d1y, d2x, d2y;
    getperpdirection(Line1,d1x,d1y); getperpdirection(Line2,d2x,d2y);
    if ((fabs(d1x-d2x)<tol) && (fabs(d1y-d2y)<tol))
      return true;
    return false;
  }

  // Get the direction of the line. This is basicaly just the normalised line such that sqrt(a^{2}+b^{2})=1.0. Note: This will be the perpendicular direction.
  template <class Line2DT, class T>
  static void getperpdirection(const Line2DT &Line, T &dx, T &dy) {
    T norm;
    norm=hypot(Line[0], Line[1]);
    if (norm!=0.0)
    { dx=Line[0]/norm; dy=Line[1]/norm; }
    else
    { dx=0.0; dy=0.0; }
    return;
  }

  /// To intersect with another line. Returns false if they don't intersect, or true if they do intersect. If they do intersect, the supplied point will be alterd to contain the intersecion point. If they don't intersect it won't be touched.
  template <class Point2DT, class Line2DT1, class Line2DT2, class T>
  static bool Intersect(Point2DT &IPoint, const Line2DT1 &Line1, const Line2DT2 &Line2, const T tol) {
    if (Parallel(Line1, Line2, tol))
      return false;
    IPoint[0]=(Line1[1]*Line2[2]-Line1[2]*Line2[1])/(Line2[1]*Line1[0]-Line1[1]*Line2[0]);
    IPoint[1]=(Line1[0]*Line2[2]-Line1[2]*Line2[0])/(Line1[1]*Line2[0]-Line1[0]*Line2[1]);
    return true;
  }

  /// Intersect a given line with a line segment specified by two points. Returns false if it doesn't intersect, otherwise IntersectP is returned containing the intersection point.
  template <class Line2DT, class Point2DT1, class Point2DT2, class Point2DT3>
  static bool IntersectWithLineSeg(const Point2DT1 &SegStart, const Point2DT2 &SegEnd, const Line2DT &Line, Point2DT3 &IntersectP) {
    // Get line between SegStart->SegEnd
    Geometry::Line2D<double> SegLine(
      SegStart.y()-SegEnd.y(),
      SegEnd.x()-SegStart.x(),
      (SegStart.x()-SegEnd.x())*SegStart.y()+(SegEnd.y()-SegStart.y())*SegStart.x());

    Geometry::Point2D<double> IPoint;
    if (!Intersect(IPoint, SegLine, Line, 1e-12))
      return false;

//    if (IPoint[0]>=MIN(SegStart[0], SegEnd[0]) && IPoint[0]<=MAX(SegStart[0], SegEnd[0]) && IPoint[1]>=MIN(SegStart[1], SegEnd[1]) && IPoint[1]<=MAX(SegStart[1], SegEnd[1]))
    if (IPoint[0]-MIN(SegStart[0], SegEnd[0])>-1e-6 && IPoint[0]-MAX(SegStart[0], SegEnd[0])<1e-6 && IPoint[1]-MIN(SegStart[1], SegEnd[1])>-1e-6 && IPoint[1]-MAX(SegStart[1], SegEnd[1])<1e-6)
    {
      IntersectP[0]=IPoint[0];
      IntersectP[1]=IPoint[1];
      return true;
    }
    return false;
  }

  // Intersect the given line (Line) with the box defined by the four points (TopLeft, TopRight, BottomRight, BottomLeft). Returns true if it does intersect and false if it doesn't. If it does intersect, the start and end points of the intersected line segment are returned in Start and End
  template <class Point2DT1, class Point2DT2, class Point2DT3, class Point2DT4, class Point2DT5, class Point2DT6, class Line2DT>
  static bool IntersectWithBox(const Point2DT1 &TopLeft, const Point2DT2 &TopRight, const Point2DT3 &BottomLeft, const Point2DT4 &BottomRight, const Line2DT &Line, Point2DT5 &Start, Point2DT6 &End) {
    // So, intersect with all the line segments
    Geometry::Point2D<double> IPoints[4];
    unsigned int i=0;

    if (IntersectWithLineSeg(TopLeft, TopRight, Line, IPoints[i]))
      i++;
    if (IntersectWithLineSeg(TopRight, BottomRight, Line, IPoints[i]))
      i++;
    if (IntersectWithLineSeg(BottomRight, BottomLeft, Line, IPoints[i]))
      i++;
    if (IntersectWithLineSeg(BottomLeft, TopLeft, Line, IPoints[i]))
      i++;

    if (i==0)
      return false;

    // Presumably intersects almost exactly with a corner (hopefuly!)
    // return start and end as the same
    if (i==1) {
      Start=IPoints[0]; End=IPoints[0];
    }

    if (i==2) {
      Start=IPoints[0]; End=IPoints[1];
    }

    if (i>=3) {
      Start=IPoints[0];
      double dist=0.0;
      for (unsigned int c=1;c<i;++c)
	if (distance_points(IPoints[c], Start)>dist)
	{ End=IPoints[c]; dist=distance_points(IPoints[c], Start); }
    }

    return true;
  }











  // --------------------------------------------
  // Functions operating in the unrectified image
  // --------------------------------------------

  // This adds a value has been added to an angle representing a polar co-ordinate. This function checks to see if the angle has looped around from -180 to 180 or vice-versa and corrects it appropriately.
  double GeneralPlanarRectify::addtoangle(double angle, const double addon) const {
    assert(addon<=2*M_PI);
    angle+=addon;
    while (angle<-M_PI)
      angle+=2.0*M_PI;
    while (angle>M_PI)
      angle-=2.0*M_PI;
   return angle;
  }

  double GeneralPlanarRectify::GetEpiDist(const Homg2DPoint<double> &e, const unsigned int xs, const unsigned int ys) const
  {
    // Check for a numerically infinite epipole
    if (e.w()>DBL_EPSILON*(MIN(fabs(e.x()), fabs(e.y()))))
      return DBL_MAX;
    Point2D<double> en;
    en.x()=e.x()/e.w(); en.y()=e.y()/e.w();
    return distance_points(Point2D<double>(xs/2.0, ys/2.0), en);
  }

  // Gets an oriented epipolar line from a point, that is to say restricted up to a positive scale factor. This function works fine even if the epipole is infinite.
  Line2D<double> GeneralPlanarRectify::GetOrientELine(const Homg2DPoint<double> &Point, const Homg2DPoint<double> &e) const {
    Line2D<double> Line;
    CrossProduct(Point, e, Line);
    // Now, give the line a sign based on quadrant the point is in in relation to the
    // epipole
    short int sign=1;
    if ((Point.x()*e.w()<e.x() && Line.b()>0.0) ||
	(Point.x()*e.w()>e.x() && Line.b()<0.0))
      sign=-1;
    if ((Point.y()*e.w()<e.y() && Line.a()<0.0) ||
	(Point.y()*e.w()>e.y() && Line.a()>0.0))
      sign=-1;
    Line.a()=Line.a()*sign;
    Line.b()=Line.b()*sign;
    Line.c()=Line.c()*sign;
    return Line;
  }

  Geometry::Point2D<double> GeneralPlanarRectify::TransferPointUsingH(const MatVec::Matrix<double> &H, const double &px, const double &py) const {
    double tz=H(2,0)*px+H(2,1)*py+H(2,2);
    return Geometry::Point2D<double>((H(0,0)*px+H(0,1)*py+H(0,2))/tz,
				     (H(1,0)*px+H(1,1)*py+H(1,2))/tz);
  }

  // Transfer a line using a homography
  Geometry::Line2D<double> GeneralPlanarRectify::TransferLineUsingH(const MatVec::Matrix<double> &H, const Geometry::Line2D<double> &Line) const {
    return Geometry::Line2D<double>(Line[0]*H(0,0)+Line[1]*H(1,0)+Line[2]*H(2,0), Line[0]*H(0,1)+Line[1]*H(1,1)+Line[2]*H(2,1), Line[0]*H(0,2)+Line[1]*H(1,2)+Line[2]*H(2,2));
  }


















  // -------------------------------------------
  // Functions operating in the nonlinear image.
  // -------------------------------------------

  double GeneralPlanarRectify::GetPointAngle(const Point2D<double> &Point) const {
    // Get the relevant half epipolar line and do it from that instead
    // That way will still work even if epipole is infinite
    Line2D<double> Line=GetOrientELine(Homg2DPoint<double>(Point.x(), Point.y(), 1.0), epipole);
    return atan2(Line.a(), Line.b());
  }

  double GeneralPlanarRectify::GetDist(const Point2D<double> &Point) const {
    return distance_points(Point, epipole);
  }

  Line2D<double> GeneralPlanarRectify::GetEpipolarLineFromAng(const double angle) const {
    return Line2D<double>(sin(angle)*epipole.w(),
			  cos(angle)*epipole.w(),
			  -epipole.x()*sin(angle)-epipole.y()*cos(angle));
  }

  Point2D<double> GeneralPlanarRectify::GetPointFromDistAndAng(const double dist, const double angle) const {
    // Now, get point at distance dist from the epipole along the line Line
    const Line2D<double> Line=GetEpipolarLineFromAng(angle);
    const double norm=hypot(Line[1],Line[0]);
    // Get normalised line gradient into dx,dy
    const double dx=Line[1]/norm, dy=-Line[0]/norm;
    return Point2D<double>(epipole.x()+dx*dist, epipole.y()+dy*dist);
  }

  Point2D<double> GeneralPlanarRectify::RectPointNLin(const Point2D<double> &PointIn) const {
    // So, first need to find the relevant angle
    const double theta=GetPointAngle(PointIn);

    // Now, find the relevant pair of surrounding angles in the map
    RectifyTabT::const_iterator first=RectTab.begin(), second=first; ++second;
    // Need the final fabs(addtoanglestuff)>M_PI to prevent problems at the wrap around from M_PI to -M_PI.
    bool iterate=true;
    while (second!=RectTab.end() && iterate) {
      if ((SGN(addtoangle(second->first, -theta))==SGN(addtoangle(first->first, -theta)) || fabs(addtoangle(second->first, -theta)-addtoangle(first->first, -theta))>M_PI)) {
	++second; ++first;
      }
      else
	iterate=false;
    }
    double angle1=MIN(first->first, second->first),
      angle2=MAX(first->first, second->first);
    double y1=first->second;

    double x=GetDist(PointIn);

    // use angle of line (theta) to interpolate an exact y value
    return Point2D<double>(x, ((theta-angle1)/(angle2-angle1))+y1);
  }

  Point2D<double> GeneralPlanarRectify::UnRectPointNLin(const Point2D<double> &PointIn) const {
    double dist=PointIn.x(), angle;

    // NonLinear Y scaling
    // So, use the y co-ordinate to look up the relevant epipolar line
    const double ryf=floor(PointIn.y());
    const int ry=(int)ryf;

    // If y is effectively integer, just read out of the table
    if (fabs(PointIn.y()-ryf)<1e-6)
      angle=UnRectTab[ry];
    else
    {
      // y is sub-pixel so need to interpolate
      double theta1=UnRectTab[ry], theta2=UnRectTab[ry+1];
      // So, to get actual line, interpolate between angles
      angle=((theta2-theta1)*(PointIn.y()-floor(PointIn.y())))+theta1;
    }
    return GetPointFromDistAndAng(dist, angle);
  }
















  // ------------------------------------------------
  // Major functions that actualy control the process
  // ------------------------------------------------
#ifdef DEBUG_GENERALRECTIFY
  template <class Iterator>
  double TransferError(Matrix<double> &H, Iterator MatchesBegin, Iterator MatchesEnd) {
    double totalerr=0.0;
    double npoints=0.0;
    for (Iterator iter=MatchesBegin; iter!=MatchesEnd; ++iter) {
      const double px=iter->first[0], py=iter->first[1], pw=1.0;
      double tx=H(0,0)*px+H(0,1)*py+H(0,2)*pw;
      double ty=H(1,0)*px+H(1,1)*py+H(1,2)*pw;
      double tz=H(2,0)*px+H(2,1)*py+H(2,2)*pw;
      Point2D<double> NP(tx/tz,ty/tz);
      totalerr+=distance_points(NP, iter->second);
      ++npoints;
    }
    return totalerr/npoints;
  }
#endif

  void GeneralPlanarRectify::SetUpTablesAndXBounds(void) {
    minxim1=DBL_MAX; maxxim1=-DBL_MAX; minxim2=DBL_MAX; maxxim2=-DBL_MAX;
    RectTab.clear(); UnRectTab.clear();

    Line2D<double> Line;
    Point2D<double> start,end;

    Point2D<double> TopLeft1(0,0), TopRight1(im1xs,0);
    Point2D<double> BottomRight1(im1xs,im1ys), BottomLeft1(0,im1ys);
    Point2D<double> TopLeft2(0,0), TopRight2(im1xs,0);
    Point2D<double> BottomRight2(im1xs,im1ys), BottomLeft2(0,im1ys);

    TopLeft1=TransferPointUsingH(Plane0Im1, TopLeft1.x(), TopLeft1.y());
    TopRight1=TransferPointUsingH(Plane0Im1, TopRight1.x(), TopRight1.y());
    BottomLeft1=TransferPointUsingH(Plane0Im1, BottomLeft1.x(), BottomLeft1.y());
    BottomRight1=TransferPointUsingH(Plane0Im1, BottomRight1.x(), BottomRight1.y());
    TopLeft2=TransferPointUsingH(Plane0Im2, TopLeft2.x(), TopLeft2.y());
    TopRight2=TransferPointUsingH(Plane0Im2, TopRight2.x(), TopRight2.y());
    BottomLeft2=TransferPointUsingH(Plane0Im2, BottomLeft2.x(), BottomLeft2.y());
    BottomRight2=TransferPointUsingH(Plane0Im2, BottomRight2.x(), BottomRight2.y());

    double angle=minang, angle2=minang;
    double maxang2=maxang;
    if (maxang<minang)
      maxang2+=2.0*M_PI;
    double x,x2;
    double angstep;

    while (angle2<maxang2) {
      // Handle the line at the current angle
      // Add it into the rectification tables
      RectTab.push_back(make_pair(angle, UnRectTab.size()));
      UnRectTab.push_back(angle);

      // Update max and min x for both images
      // Do image 1 first
      Line=GetEpipolarLineFromAng(angle);
      // Update max and min x for line
      IntersectWithBox(TopLeft1, TopRight1, BottomLeft1, BottomRight1,
        Line, start, end);
      x=GetDist(start); x2=GetDist(end);
      minxim1=MIN(x, minxim1); maxxim1=MAX(x, maxxim1);
      minxim1=MIN(x2, minxim1); maxxim1=MAX(x2, maxxim1);
      // Get maximum angle step for image 1 that preserves worst case pixel size
      angstep=atan(1.0/MAX(x,x2));

      // Now image 2
      // Already have the epipolar line in image 1
      // Update max and min x for line
      IntersectWithBox(TopLeft2, TopRight2, BottomLeft2, BottomRight2,
        Line, start, end);

      x=GetDist(start); x2=GetDist(end);
      minxim2=MIN(x, minxim2); maxxim2=MAX(x, maxxim2);
      minxim2=MIN(x2, minxim2); maxxim2=MAX(x2, maxxim2);

      // Make angstep minimal step
      angstep=MIN(angstep, atan(1.0/MAX(x,x2)));

      // Now, find the next epipolar line
      // Is smallest angle step of step in image 2 and image 1
      angle=addtoangle(angle, angstep);
      angle2+=angstep;
    }

    // Make min and max x same for both images so both images are the same size
    minxim1=MIN(minxim1, minxim2);
    maxxim1=MAX(maxxim1, maxxim2);
    minxim2=minxim1; maxxim2=maxxim1;

#ifdef DEBUG_GENERALRECTIFY
    cerr << "Image 1: minx=" << minxim1 << " maxx=" << maxxim1 << endl;
    cerr << "Image 2: minx=" << minxim2 << " maxx=" << maxxim2 << endl;
    cerr << "Output image y sizes = " << (unsigned int)UnRectTab.size() << endl;
#endif
  }

  void GeneralPlanarRectify::GetImAngleRange(const unsigned int xs, const unsigned int ys, double &minangim, double &maxangim, const Matrix<double> &H) const {
    Point2D<double> TopLeft(0,0), TopRight(xs,0);
    Point2D<double> BottomRight(xs,ys), BottomLeft(0,ys);

    TopLeft=TransferPointUsingH(H, TopLeft.x(), TopLeft.y());
    TopRight=TransferPointUsingH(H, TopRight.x(), TopRight.y());
    BottomLeft=TransferPointUsingH(H, BottomLeft.x(), BottomLeft.y());
    BottomRight=TransferPointUsingH(H, BottomRight.x(), BottomRight.y());

    // Check where the epipole is.
    // So, check if epipole is within transfered image by intersecting
    // line with each corner with the image box. If all lines intersect the box
    // twice then is within the image
    Line2D<double> TopL, TopR, BotL, BotR;
    CrossProduct(epipole, Homg2DPoint<double>(TopLeft.x(), TopLeft.y(), 1.0), TopL);
    CrossProduct(epipole, Homg2DPoint<double>(TopRight.x(), TopRight.y(), 1.0), TopR);
    CrossProduct(epipole, Homg2DPoint<double>(BottomLeft.x(), BottomLeft.y(), 1.0), BotL);
    CrossProduct(epipole, Homg2DPoint<double>(BottomRight.x(), BottomRight.y(), 1.0), BotR);

    Point2D<double> start1,end1;
    Point2D<double> start2,end2;
    Point2D<double> start3,end3;
    Point2D<double> start4,end4;
    if (IntersectWithBox(TopLeft, TopRight, BottomLeft, BottomRight,
				    TopL, start1, end1) &&
	IntersectWithBox(TopLeft, TopRight, BottomLeft, BottomRight,
				    TopR, start2, end2) &&
	IntersectWithBox(TopLeft, TopRight, BottomLeft, BottomRight,
				    BotL, start3, end3) &&
	IntersectWithBox(TopLeft, TopRight, BottomLeft, BottomRight,
				    BotR, start4, end4)) {
      // If all lines intersected with the box twice then the epipole must be within
      // the image. Need to be careful of rounding errors, so don't check exact
      // equivalence
      double THRESH=1e-4;
      if (distance_points(start1,end1)>THRESH &&
	  distance_points(start2,end2)>THRESH &&
	  distance_points(start3,end3)>THRESH &&
	  distance_points(start4,end4)>THRESH) {
	minangim=-M_PI; maxangim=M_PI; return;
      }
    }

    double addon=-GetPointAngle(TopLeft), angle;
    minangim=0.0; maxangim=0.0;
    angle=GetPointAngle(TopRight); angle=addtoangle(angle, addon);
    minangim=MIN(minangim, angle); maxangim=MAX(maxangim, angle);
    angle=GetPointAngle(BottomRight); angle=addtoangle(angle, addon);
    minangim=MIN(minangim, angle); maxangim=MAX(maxangim, angle);
    angle=GetPointAngle(BottomLeft); angle=addtoangle(angle, addon);
    minangim=MIN(minangim, angle); maxangim=MAX(maxangim, angle);

    minangim=addtoangle(minangim, -addon);
    maxangim=addtoangle(maxangim, -addon);
  }


























  //-------------------------------------
  // Public functions for use by the user
  //-------------------------------------

  GeneralPlanarRectify::GeneralPlanarRectify(const unsigned int im1xs_, const unsigned int im1ys_,
					     const unsigned int im2xs_, const unsigned int im2ys_,
					     const FMatrix<double> &FM_,
					     const vector<pair<Point2D<double>, Point2D<double> > > &Matches_)
    : im1xs(im1xs_), im1ys(im1ys_),
      im2xs(im2xs_), im2ys(im2ys_), FM(FM_) {

    Homg2DPoint<double> e1, e2;
    GetEpipoles(FM, e1, e2);

#ifdef DEBUG_GENERALRECTIFY
    cerr << "Image sizes: " << im1xs << "x" << im1ys << endl;
    cerr << "Image 1 epipole: (" << e1.x()/e1.w() << "," << e1.y()/e1.w() << ")";
    cerr << " Image 2 epipole: (" << e2.x()/e2.w() << "," << e2.y()/e2.w() << ")" << endl;
#endif

    typedef vector<pair<Point2D<double>, Point2D<double> > > MatchVec2Im;
    // Make all the matches perfect
    MatchVec2Im Matches, OutlierFree;
    for (unsigned int i=0;i<Matches_.size();++i)
//      Matches.push_back(FM.get_nearest_perfect_match(Matches_[i], 1e-7));
      Matches.push_back(Matches_[i]);

    // ---------------------------------------------
    // Now, get the consistant homography (robustly)
    // ---------------------------------------------
    Matrix<double> H(3,3), HInv(3,3);

    // Yeuch - what WAS I thinking when I wrote the calcplanarhomog routines.
    // I did them the wrong way around :o(
    HomogConsistF::RANSAC(Matches, OutlierFree, HInv, FM);
    MatrixOps::inverse(HInv, H);

#ifdef DEBUG_GENERALRECTIFY
    cerr << "Got Full Homography, transfer error=" << TransferError(H, OutlierFree.begin(), OutlierFree.end()) << " outliers=" << (unsigned int)(Matches.size()-OutlierFree.size()) << endl;
/*
    Geometry:: Point2D<double> e1n(e1.x()/e1.w(), e1.y()/e1.w()), e2n(e2.x()/e2.w(), e2.y()/e2.w());
    MultiViewGeom::EuclMatch2Im ematch(e2n,e1n);
    cerr << TransferError(H, &ematch, &ematch+1) << endl;*/
#endif

    Matrix<double> id(3,3);
    for (size_t i=0;i<3;++i)
      for (size_t j=0;j<3;++j) {
	if (i==j)
	  id(i,j)=1.0;
	else
	  id(i,j)=0.0;
      }
    // Select image with epipole closest to the image as the nonlinear image
    if (GetEpiDist(e1, im1xs, im1ys)<GetEpiDist(e2, im2xs, im2ys)) {
      epipole=e1;
      Plane0Im1=id; Plane0InvIm1=Plane0Im1;
      Plane0Im2=H; Plane0InvIm2=HInv;
#ifdef DEBUG_GENERALRECTIFY
      cerr << "Doing nonlinear transform in image 1" << endl;
#endif
    }
    else
    {
      epipole=e2;
      Plane0Im1=HInv; Plane0InvIm1=H;
      Plane0Im2=id; Plane0InvIm2=Plane0Im2;
#ifdef DEBUG_GENERALRECTIFY
      cerr << "Doing nonlinear transform in image 2" << endl;
#endif
    }
    epipole.x()=epipole.x()/epipole.w();
    epipole.y()=epipole.y()/epipole.w();
    epipole.w()=1.0;

    // ------------------------------------------
    // Get the common angle range for both images
    // ------------------------------------------
    // First get extreme angles from both images
    double minangim1, maxangim1, minangim2, maxangim2;
    GetImAngleRange(im1xs, im1ys, minangim1, maxangim1, Plane0Im1);
    GetImAngleRange(im2xs, im2ys, minangim2, maxangim2, Plane0Im2);
#ifdef DEBUG_GENERALRECTIFY
    cerr << "  Image 1 min angle=" << minangim1 << " max angle=" << maxangim1 << endl;
    cerr << "  Image 2 min angle=" << minangim2 << " max angle=" << maxangim2 << endl;
#endif

    // If the epipole is in either image, then have anomalous case
    if (maxangim1-minangim1==2.0*M_PI || maxangim2-minangim2==2.0*M_PI)
    {
      // Set to range of image without epipole in it
      // If both within then will simply be full range anyway
      if (maxangim1-minangim1==2.0*M_PI)
      { minang=minangim2; maxang=maxangim2; }
      else
      { minang=minangim1; maxang=maxangim1; }
    }
    else
    {
      // Otherwise need to take second minimum and second maximum angle
      if (fabs(minangim1-minangim2)>M_PI)
	minang=MIN(minangim1, minangim2);
      else
	minang=MAX(minangim1, minangim2);

      if (fabs(maxangim1-maxangim2)>M_PI)
	maxang=MAX(maxangim1, maxangim2);
      else
	maxang=MIN(maxangim1, maxangim2);
    }
#ifdef DEBUG_GENERALRECTIFY
    cerr << "  Common range, min angle=" << minang << " max angle=" << maxang << endl;
#endif

    SetUpTablesAndXBounds();

    // Now, check to see if need to flip the images about y axis
    // Do flips by flipping the angle range (add 180 degrees), and
    // using -ve distances
    // Purely aesthetic and totaly unimportant
    flipx=false;
    if (SGN(e1.x())*SGN(e1.w())==1)
      flipx=true;
#ifdef DEBUG_GENERALRECTIFY
    cerr << "flipx: ";
    if (flipx)
      cerr << "true" << endl;
    else
      cerr << "false" << endl;
#endif

  }

  GeneralPlanarRectify::~GeneralPlanarRectify() {
  }

  void GeneralPlanarRectify::RectifyPointIm1(double &x, double &y) const {
    // First, convert the point to the nonlinear image
    Point2D<double> PointIn(TransferPointUsingH(Plane0Im1, x, y));

    PointIn=RectPointNLin(PointIn);
    // And account for image bounds and flipping
    PointIn.x()-=minxim1;
    if (flipx)
      PointIn.x()=(maxxim1-minxim1)-PointIn.x();
    x=PointIn.x(); y=PointIn.y();
  }
  void GeneralPlanarRectify::RectifyPointIm2(double &x, double &y) const {
    // First, convert the point to the nonlinear image
    Point2D<double> PointIn(TransferPointUsingH(Plane0Im2, x, y));

    PointIn=RectPointNLin(PointIn);
    // And account for image bounds and flipping
    PointIn.x()-=minxim2;
    if (flipx)
      PointIn.x()=(maxxim2-minxim2)-PointIn.x();
    x=PointIn.x(); y=PointIn.y();
  }

  void GeneralPlanarRectify::UnRectifyPointIm1(double &x, double &y) const {
    // First, remove any bounding and flipping effects
    if (flipx)
      x=maxxim1-x;
    else
      x+=minxim1;

    // Undo nonlinear rectification
    Point2D<double> PointIn(UnRectPointNLin(Point2D<double>(x,y)));

    // And transfer to image one
    PointIn=TransferPointUsingH(Plane0InvIm1, PointIn.x(), PointIn.y());    

    x=PointIn.x(); y=PointIn.y();
  }
  void GeneralPlanarRectify::UnRectifyPointIm2(double &x, double &y) const {
    // First, remove any bounding and flipping effects
    if (flipx)
      x=maxxim2-x;
    else
      x+=minxim2;

    // Undo nonlinear rectification
    Point2D<double> PointIn(UnRectPointNLin(Point2D<double>(x,y)));

    // And transfer to image one
    PointIn=TransferPointUsingH(Plane0InvIm2, PointIn.x(), PointIn.y());
    x=PointIn.x(); y=PointIn.y();
  }

  // This function unrectifies a set of points in image 1 all on the same epipolar
  // line
  void GeneralPlanarRectify::UnRectifyPointsIm1(vector<double> xvals,
						const double y,
						vector<Point2D<double> > &RectPoints) const {

    double angle;
    // So, get the angle of the unrectified line first
    // NonLinear Y scaling
    // So, use the y co-ordinate to look up the relevant epipolar line
    const double ryf=floor(((double)y));
    const int ry=(int)ryf;

    // If y is effectively integer, just read out of the table
    if (fabs(y-ryf)<1e-6)
      angle=UnRectTab[ry];
    else
    {
      // y is sub-pixel so need to interpolate
      double theta1=UnRectTab[ry], theta2=UnRectTab[ry+1];
      // So, to get actual line, interpolate between angles
      angle=((theta2-theta1)*(y-floor(y)))+theta1;
    }
  
    const Line2D<double> ELine(GetEpipolarLineFromAng(angle));
    const double norm=hypot(ELine[1],ELine[0]);
    // Get normalised line gradient into dx,dy
    const double dx=ELine[1]/norm, dy=-ELine[0]/norm;

    double dist;
    vector<double>::iterator BeginX=xvals.begin(), EndX=xvals.end();
    for (; BeginX!=EndX; ++BeginX) {
      // So, get dist for the given point
      if (flipx)
	dist=maxxim1-*BeginX;
      else
	dist=*BeginX+minxim1;
      // And find point at distance dist from the epipole on the given line
      Point2D<double> Im1P(epipole.x()+dx*dist, epipole.y()+dy*dist);
      RectPoints.push_back(TransferPointUsingH(Plane0InvIm1, Im1P.x(), Im1P.y()));
    }
  }

  // This function unrectifies a set of points in image 2 all on the same epipolar
  // line
  void GeneralPlanarRectify::UnRectifyPointsIm2(vector<double> xvals,
					      const double y,
					      vector<Point2D<double> > &RectPoints) const {
    double angle;
    // So, get the angle of the unrectified line first
    // NonLinear Y scaling
    // So, use the y co-ordinate to look up the relevant epipolar line
    const double ryf=floor(y);
    const int ry=(int)ryf;

    // If y is effectively integer, just read out of the table
    if (fabs(y-ryf)<1e-6)
      angle=UnRectTab[ry];
    else
    {
      // y is sub-pixel so need to interpolate
      double theta1=UnRectTab[ry], theta2=UnRectTab[ry+1];
      // So, to get actual line, interpolate between angles
      angle=((theta2-theta1)*(y-floor(y)))+theta1;
    }
  
    const Line2D<double> ELine(GetEpipolarLineFromAng(angle));
    const double norm=hypot(ELine[1],ELine[0]);
    // Get normalised line gradient into dx,dy
    const double dx=ELine[1]/norm, dy=-ELine[0]/norm;

    vector<double>::iterator BeginX=xvals.begin(), EndX=xvals.end();
    for (; BeginX!=EndX; ++BeginX) {
      // So, get dist for the given point
      double dist;
      if (flipx)
	dist=maxxim2-*BeginX;
      else
	dist=*BeginX+minxim2;
      // And find point at distance dist from the epipole on the given line
      Point2D<double> Im2P(epipole.x()+dx*dist, epipole.y()+dy*dist);
      RectPoints.push_back(TransferPointUsingH(Plane0InvIm2, Im2P.x(), Im2P.y()));
    }
  }

/*
  std::ostream& GeneralPlanarRectify::SaveToStream(std::ostream &s) const {
    s << "General Planar Rectification:" << endl;
    s << epipole << FM << Plane0Im1 << Plane0InvIm1 << Plane0Im2 << Plane0InvIm2;
    s << im1xs << " " << im1ys << " " << im2xs << " " << im2ys << " ";
    s << minang << " " << maxang << " ";
    s << minxim1 << " " << maxxim1 << " " << minxim2 << " " << maxxim2 << " " << flipx << endl;
    return s;
  }

  std::istream& GeneralPlanarRectify::LoadFromStream(std::istream &s) {
    string LookFor("General Planar Rectification:");
    std::string::size_type LForSize=LookFor.size();
    char Temp[LForSize];
    int off=s.tellg();
    s.read(Temp, LForSize);
    if (s.eof())
      s.seekg(off, std::ios::beg);

    if (LookFor.compare(Temp)!=0) {
      s.seekg(-((int)LForSize), std::ios::cur);
      s.setstate(ios_base::badbit);
      return s;
    }
    s >> epipole >> FM >> Plane0Im1 >> Plane0InvIm1 >> Plane0Im2 >> Plane0InvIm2 >> im1xs >> im1ys >> im2xs >> im2ys >> minang >> maxang >> minxim1 >> maxxim1 >> minxim2 >> maxxim2 >> flipx; s.ignore(1);
    SetUpTablesAndXBounds();
    return s;
  }
*/

} // namespace Rectify
