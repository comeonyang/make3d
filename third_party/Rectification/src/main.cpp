#include <iostream>
#include <fstream>
#include <time.h> // For speed checks
#include <float.h>
#include <vector>
#include <iomanip> // For machine precision

#include "main.h"
#include "PNM.h"
#include "Rectify/General2Im.h"

// Keep the fortran 2 C libraries happy. (Give the fortran main function).
extern "C" {
  void MAIN__() {}
}

/*********************/
/* The main function */
/*********************/
int main(int argc, char **argv)
{
  // Default exit status to ok.
  int ExitStatus=0;

  if (argc!=7) {
    std::cerr << "syntax: " << argv[0] << " image1 image2 fmatrix matches outimage1 outimage2" << std::endl;
    std::cerr << std::endl << "where:" << std::endl;
    std::cerr << "    image1, image2 Are the two images to be rectified" << std::endl;
    std::cerr << "    fmatrix is the fundamental matrix" << std::endl;
    std::cerr << "    outimage1, outimage2 are filenames for rectified images" << std::endl;
    exit(1);
  }

  // Load images
  Image::Image<unsigned int> image1, image2, rectim1, rectim2;

  PNMLoadSave::LoadImage(argv[1], image1);
  PNMLoadSave::LoadImage(argv[2], image2);

  // Load the fundamental matrix
  MultiViewGeom::FMatrix<double> FM;
  std::ifstream FMfile(argv[3], std::ios::in | std::ios::binary);
  if (!FMfile) {
    std::cerr << "couldn't open fundamental matrix file " << argv[3] << std::endl;
    exit(1);
  }
  FMfile >> std::setprecision(DBL_DIG);
  for (int row=0;row<3;++row) {
    for (int col=0;col<3;++col) {
      FMfile >> FM(row,col);
    }
  }

  std::vector<std::pair<Geometry::Point2D<double>, Geometry::Point2D<double> > > Matches;
  std::ifstream Matchesfile(argv[4], std::ios::in | std::ios::binary);
  if (!Matchesfile) {
    std::cerr << "couldn't open file containing matches " << argv[4] << std::endl;
    exit(1);
  }
  Matchesfile >> std::setprecision(DBL_DIG);
  while (!Matchesfile.eof()) {
    double x,y,x2,y2;
    Matchesfile >> x >> y >> x2 >> y2;
    if (!Matchesfile.eof())
      Matches.push_back(std::make_pair(Geometry::Point2D<double>(x,y), Geometry::Point2D<double>(x2,y2)));
  }

  Rectify::GeneralPlanarRectify rectify(image1.xsize(), image1.ysize(), image2.xsize(), image2.ysize(), FM, Matches);
  rectify.resampleIms(image1, image2, rectim1, rectim2, 0xFF000000);

  PNMLoadSave::SaveImage(argv[5], rectim1);
  PNMLoadSave::SaveImage(argv[6], rectim2);

  return ExitStatus;
}
