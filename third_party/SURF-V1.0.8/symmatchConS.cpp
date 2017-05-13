/*
 * Speeded-Up Robust Features (SURF)
 * http://people.ee.ethz.ch/~surf
 *
 * Sample application for feature matching using nearest-neighbor
 * ratio method.
 *
 * BUILD USING "make match.ln".
 *
 * Author: Andreas Ess
 *
 * Copyright (2006): ETH Zurich, Switzerland
 * Katholieke Universiteit Leuven, Belgium
 * All rights reserved.
 *
 * For details, see the paper:
 * Herbert Bay,  Tinne Tuytelaars,  Luc Van Gool,
 *  "SURF: Speeded Up Robust Features"
 * Proceedings of the ninth European Conference on Computer Vision, May 2006
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for educational, research, and non-commercial
 * purposes, without fee and without a signed licensing agreement, is
 * hereby granted, provided that the above copyright notice and this
 * paragraph appear in all copies modifications, and distributions.
 *
 * Any commercial use or any redistribution of this software
 * requires a license from one of the above mentioned establishments.
 *
 * For further details, contact Andreas Ess (aess@vision.ee.ethz.ch).
 */

/**
   modified by Jeff Michels to add the requirement that matches
   are only output if the point in im2 is the best match for the 
   point in im1 AND vice versa.
**/

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "ipoint.h"


#include "image.h"
#include "imload.h"


using namespace std;
using namespace surf;

// Length of descriptor vector, ugly, global variable
unsigned int vlen;

// Calculate square distance of two vectors
double distSquare(double *v1, double *v2, int n) {
	double dsq = 0.;
	while (n--) {
		dsq += (*v1 - *v2) * (*v1 - *v2);
		v1++;
		v2++;
	}
	return dsq;
}

// Find closest interest point in a list, given one interest point
int findMatch(const Ipoint& ip1, const vector< Ipoint >& ipts) {
	double mind = 1e100, second = 1e100; //Magic number Min??
	int match = -1;

	for (unsigned i = 0; i < ipts.size(); i++) {
		// Take advantage of Laplacian to speed up matching
		if (ipts[i].laplace != ip1.laplace)
			continue;

		double d = distSquare(ipts[i].ivec, ip1.ivec, vlen);

		if (d < mind) {
			second = mind;
			mind = d;
			match = i;
		} else if (d < second) {
			second = d;
		}
	}

	if (mind < 0.7 * second) //0.8 //Magic number Min?? 
            //Answer: in SIFT paper 0.8 is used for discarding 90% of false match but discard only 5% of correct matches
		return match;

	return -1;
}

void printAllDistances(const vector< Ipoint >& ipts1, const vector< Ipoint >& ipts2) {
  for (unsigned int i = 0; i < ipts1.size(); i++) {
    for (unsigned int j = 0; j < ipts2.size(); j++) {
      double dist = distSquare(ipts1[i].ivec, ipts2[j].ivec, vlen);
      cout << i << " " << j << " " << dist << endl;
    }
  }
}

//Min add ConstrainCheck to check if the Constrain been meet
bool ConstrainCheck(const Ipoint& ipts1, const Ipoint& ipts2, const Ipoint& ConS1, const Ipoint& ConS2){

//	cout << "test ConS1" << ConS1.ivec[1];
//	cout << "test ConS2" << ConS2.ivec[1];
	
	bool Sat = false;
	// check if ipts2 satisfy constrain1 or ipts1 satisfy constrain2
	if (( ipts2.x >= ConS1.ivec[0] && ipts2.x <= ConS1.ivec[1] && ipts2.y >= ConS1.ivec[2] && ipts2.y <= ConS1.ivec[3]) || ( ipts1.x >= ConS2.ivec[0] && ipts1.x <= ConS2.ivec[1] && ipts1.y >= ConS2.ivec[2] && ipts1.y <= ConS2.ivec[3])){
		Sat = true;
//		cout << "violate satisfied";
	}else{
//		cout << "violate constrain";
	}
	
	return Sat;        
}

// essentially the same as findMatches except that it deletes matches that
// aren't symetric.
// Really time-consuming function
vector <int> symMatches(const vector< Ipoint >& ipts1, const vector< Ipoint >& ipts2, const vector< Ipoint >& ConS1, const vector< Ipoint >& ConS2) {
  vector< int > bestIndex1(ipts1.size());
  vector< int > bestIndex2(ipts2.size());
  vector< double > bestScore1(ipts1.size());
  vector< double > bestScore2(ipts2.size());
  vector< double > secondScore1(ipts1.size());
  vector< double > secondScore2(ipts2.size());

  // Initialize all of the vectors
  for (unsigned int i = 0; i < ipts1.size(); i++) {
    bestIndex1[i] = -1;
    bestScore1[i]= 1e100;
    secondScore1[i] = 1e100;
  }
  for (unsigned int j = 0; j < ipts2.size(); j++) {
    bestIndex2[j] = -1;
    bestScore2[j]= 1e100;
    secondScore2[j] = 1e100;
  }

  // Main loop - per point in im1 do a linear scan of points in im2
  int total = 0;
  for (unsigned int i = 0; i < ipts1.size(); i++) {
    for (unsigned int j = 0; j < ipts2.size(); j++) {
       if (ipts1[i].laplace != ipts2[j].laplace)
                        continue;
	//Min add ConstrainCheck();
	if ( !ConstrainCheck(ipts1[i], ipts2[j], ConS1[i], ConS2[j]))
                        continue;

	//record best and second best scores and best index for each
      //point in each image
      double dist = distSquare(ipts1[i].ivec, ipts2[j].ivec, vlen);
      if (dist < bestScore1[i]) {
	bestIndex1[i] = j;
	secondScore1[i] = bestScore1[i];
	bestScore1[i] = dist;
      } else if (dist < secondScore1[i]) {
	secondScore1[i] = dist;
      }
      if (dist < bestScore2[j]) {
	bestIndex2[j] = i;
	secondScore2[j] = bestScore2[j];
	bestScore2[j] = dist;
      } else if (dist < secondScore2[j]) {
	secondScore2[j] = dist;
      }
    }
  }
  //cout << "Matched " << total << " points." << endl;
  //return bestIndex1; //for threshold

  /*******************************
     //killed same img matches - should be covered by match symetry anyway
  
  //update second scores with best match from same image (for image 1)
  for (int i = 0; i < ipts1.size(); i++) {
    for (int j = i+1; j < ipts1.size(); j++) { 
      double dist = distSquare(ipts1[i].ivec, ipts1[j].ivec, vlen);
      if (dist < secondScore1[i]){// && dist < 0.3) {
	//cout << "im1 point " << i << " best: " << bestScore1[i] << " second: " << secondScore1[i] << 
	//  "im1 point " << j << ": " << dist << endl;
	secondScore1[i] = dist;
      }
    }
  }

  //update second scores with best match from same image (for image 2)
  for (int i = 0; i < ipts2.size(); i++) {
    for (int j = i+1; j < ipts2.size(); j++) { 
      double dist = distSquare(ipts2[i].ivec, ipts2[j].ivec, vlen);
      if (dist < secondScore2[j]){// && dist < 0.3){
	//cout << "im2 point " << i << " best: " << bestScore2[i] << " second: " << secondScore2[i] << 
	//  "im2 point " << j << ": " << dist << endl;
	secondScore2[i] = dist;
      }
    }
  }
  ******************************/

  //if the best score is no longer the best or isn't good enough relative
  // to the second, kill the index

  //image 1 //Min change 0.85 to 0.7
  for (unsigned int i = 0; i < ipts1.size(); i++) {
    if (secondScore1[i] * 0.7 < bestScore1[i]) {
      bestIndex1[i] = -1;
    }
  }

  //image 2 //Min change 0.85 to 0.7
  for (unsigned int i = 0; i < ipts2.size(); i++) {
    if (secondScore2[i] * 0.7 < bestScore2[i]) {
      bestIndex2[i] = -1;
    }
  }

  //if the best matches aren't mutual, kill the match
  for (unsigned int i = 0; i < ipts1.size(); i++) {
    int index = bestIndex1[i];
    if (index != -1) {
      if (bestIndex2[index] != i) {
	//cout << "Non symetric match.  im1(" << i << ")->" << index << " but im2(" << index <<
	//  ")->" << bestIndex2[index] << endl;
	bestIndex1[i] = -1;
      } else {
	cout << " Matched feature " << i << " in image 1 with feature "
	     << index << " in image 2." << endl;
	total++;
      }
    }
  }

  cout << "Matched " << total << " points." << endl;
  return bestIndex1;
}

// Find all possible matches between two images 
// need to change to reduce the search space by knowing GPS information Min
vector< int > findMatches(const vector< Ipoint >& ipts1, const vector< Ipoint >& ipts2) {
	vector< int > matches(ipts1.size());
	int c = 0;
	for (unsigned i = 0; i < ipts1.size(); i++) {
		int match = findMatch(ipts1[i], ipts2);
		matches[i] = match;
		if (match != -1) {
			cout << " Matched feature " << i << " in image 1 with feature "
				<< match << " in image 2." << endl;
			c++;
		}
	}
	cout << " --> Matched " << c << " features of " << ipts1.size() << " in image 1." << endl;
	return matches;
}

// Load the interest points from a regular ASCII file
void loadIpoints(string sFileName, vector< Ipoint >& ipts) {
	ifstream ipfile(sFileName.c_str());
	if( !ipfile ) {
		cerr << "ERROR in loadIpoints(): "
			<< "Couldn't open file '" << sFileName.c_str() << "'!" << endl;
		return;
	}

	// Load the file header
	unsigned count;
	ipfile >> vlen >> count;

	// create a new interest point vector
	ipts.clear();
	ipts.resize(count);
	//cout << "The number of features" << count << endl;

	// Load the interest points in Mikolajczyk's format
	for (unsigned n = 0; n < count; n++) {
		// circular regions with diameter 5 x scale
		float x, y, a, b, c;

		// Read in region data, though not needed for actual matching
		ipfile >> x >> y >> a >> b >> c;

		float det = sqrt((a-c)*(a-c) + 4.0*b*b);
		float e1 = 0.5*(a+c + det);
		float e2 = 0.5*(a+c - det);
		float l1 = (1.0/sqrt(e1));
		float l2 = (1.0/sqrt(e2));
		float sc = sqrt( l1*l2 );

		ipts[n].x = x;
		ipts[n].y = y;
		ipts[n].scale = sc/2.5;

		// Read in Laplacian
		ipfile >> ipts[n].laplace;

		// SURF makes Laplacian part of descriptor, so skip it..
		ipts[n].ivec = new double[vlen - 1];
		for (unsigned j = 0; j < vlen - 1; j++)
			ipfile >> ipts[n].ivec[j];
	}
}

void drawLine(Image *im, int x1, int y1, int x2, int y2) {
	if ((x1 < 0 && x2 < 0) || (y1 < 0 && y2 < 0) ||
	    (x1 >= im->getWidth() && x2 >= im->getWidth()) ||
            (y1 >= im->getHeight() && y2 >= im->getHeight()))
		return;

	bool steep = std::abs(y2 - y1) > std::abs(x2 - x1);
	if (steep) {
		int t;
		t = x1;
		x1 = y1;
		y1 = t;
		t = y2;
		y2 = x2;
		x2 = t;
	}

	if (x1 > x2) {
		// Swap points
		int t;
		t = x1;
		x1 = x2;
		x2 = t;
		t = y1;
		y1 = y2;
		y2 = t;
	}

	int deltax = x2 - x1;
	int deltay = std::abs(y2 - y1);

	int error = 0;
	int y = y1;
	int ystep = y1 < y2 ? 1 : -1;

	for (int x = x1; x < x2; x++) {
		if (steep) {
			if (x >= 0 && y >= 0 && y < im->getWidth() && x < im->getHeight())
				im->setPix(y, x, 1);
		} else {
			if (x >= 0 && y >= 0 && x < im->getWidth() && y < im->getHeight())
				im->setPix(x, y, 1);

		}
		error += deltay;

		if (2 * error > deltax) {
			y += ystep;
			error -= deltax;
		}
	}
}

void drawCross(Image *im, int x, int y, int s = 5) {
	for (int x1 = x - s; x1 <= x + s; x1++)
		im->setPix(x1, y, 1);
	for (int y1 = y - s; y1 <= y + s; y1++)
		im->setPix(x, y1, 1);
}

int main(int argc, char **argv) {
	Image *im1, *im2;

	ImLoad ImageLoader;
	vector< Ipoint > ipts1, ipts2, ConS1, ConS2; //Min 
	bool drawc = false;
	bool usesym = false;
	bool printAllDist = false;

	char ofname[100];

	im1 = im2 = NULL;
	ofname[0] = 0;

	// Read the arguments
	int arg = 0;
	while (++arg < argc) { 
		if (! strcmp(argv[arg], "-k1"))
			loadIpoints(argv[++arg], ipts1);
		if (! strcmp(argv[arg], "-k2"))
			loadIpoints(argv[++arg], ipts2);
		//Min -------------------------------------
		if (! strcmp(argv[arg], "-S1")) 
			loadIpoints(argv[++arg], ConS1);
		if (! strcmp(argv[arg], "-S2"))
			loadIpoints(argv[++arg], ConS2);
		//Min -------------------------------------

		if (! strcmp(argv[arg], "-im1"))
			im1 = ImageLoader.readImage(argv[++arg]); 
		if (! strcmp(argv[arg], "-im2"))
			im2 = ImageLoader.readImage(argv[++arg]); 
		if (! strcmp(argv[arg], "-o"))
			strcpy(ofname, argv[++arg]);
		if (! strcmp(argv[arg], "-c"))
			drawc = true;
		if (! strcmp(argv[arg], "-s"))
			usesym = true;
		if (! strcmp(argv[arg], "-d"))
			printAllDist = true;

	}

	if (ipts1.size() == 0 || ipts2.size() == 0) {
		cout << "Usage:" << endl;
		cout << " match [-s] -k1 out1.surf -k2 out2.surf -im1 img1.pgm -im2 img2.pgm -o out.pgm" << endl << endl;
		cout << "For each feature in first descriptor file, find best in second according to "
			<< "nearest neighbor ratio strategy. Display matches in out.pgm, generated "
			<< "from img1.pgm and img2.pgm. Use -c to draw crosses at interest points." << endl;

		cout << "use -s to require that the best match be symetric."<<endl;
		cout << "use -d to print out all of the n*m distances between points instead of matching."<<endl;
		return 1;
	}

	vector< int > matches;
	  
	if (usesym) {
	  cerr << "Using symmetricConS matches, using descriptor length " << vlen << ConS1[0].ivec[3] << endl;
//	  cout << "test CONS" << ConS1[1].ivec[1] << endl;
		matches = symMatches(ipts1, ipts2, ConS1, ConS2);//Min
	} else if (printAllDist){
	  printAllDistances(ipts1, ipts2);
	} else {
	  matches = findMatches(ipts1, ipts2);
	}


	if (im1 != NULL && im2 != NULL && ofname[0] != 0) {
		Image res(max(im1->getWidth(), im2->getWidth()), im1->getHeight() + im2->getHeight());
		for (int x = 0; x < im1->getWidth(); x++)
			for (int y = 0; y < im1->getHeight(); y++)
				res.setPix(x, y, im1->getPix(x, y));
		for (int x = 0; x < im2->getWidth(); x++)
			for (int y = 0; y < im2->getHeight(); y++)
				res.setPix(x, y + im1->getHeight(), im2->getPix(x, y));

		// Draw lines for matches
		for (unsigned i = 0; i < matches.size(); i++) {
			if (matches[i] != -1) {
				drawLine(&res, (int)ipts1[i].x, (int)ipts1[i].y,
					(int)ipts2[matches[i]].x, (int)(ipts2[matches[i]].y + im1->getHeight()));
			}
		}

		// Draw crosses at interest point locations
		if (drawc) {
			for (unsigned i = 0; i < ipts1.size(); i++)
				drawCross(&res, (int)ipts1[i].x, (int)ipts1[i].y);
			for (unsigned i = 0; i < ipts2.size(); i++)
				drawCross(&res, (int)ipts2[i].x, (int)ipts2[i].y + im1->getHeight());
		}

		ImageLoader.saveImage(ofname, &res);
	}

	
	return 0;
}
