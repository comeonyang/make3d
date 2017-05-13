==================== GENERAL ====================
This directory contains eucsbademo, an example of using sba for Euclidean bundle adjustment.
Refer to ICS/FORTH TR-340 for documentation on sba. eucsbademo accepts 3 file names as arguments:
They are the initial estimates for the camera motion parameters, the initial estimates for
the 3D point parameters along with the 3D points image projections and the camera intrinsic
calibration parameters. The file for the camera motion parameters has a separate line for
every camera, each line containing 7 parameters (a 4 element quaternion for rotation and a
3 element vector for translation). The file for 3D points and image projections is made up
of lines of the form 

X Y Z  NFRAMES  FRAME0 x0 y0  FRAME1 x1 y1 ...

each corresponding to a single 3D point. X Y Z are the points' Euclidean 3D coordinates,
NFRAMES the total number of images in which the points' projections are available and each
of the NFRAMES subsequent triplets FRAME x y pecifies that the 3D point in question projects
to pixel x y in image FRAME. For example, the line

100.0 200.0 300.0 3  2 270.0 114.1 4 234.2 321.7 5 173.6 425.8

specifies the 3D point (100.0, 200.0, 300.0) projecting to the 3 points (270.0, 114.1),
(234.2, 321.7) and (173.6, 425.8) in images 2, 4 and 5 respectively. Camera and 3D point
indices count from 0.

==================== FILES ====================
eucsbademo.c: main demo program
readparams.c: functions to read the initial motion and structure estimates from text files
imgproj.c:    functions to estimate the projection on a given camera of a certain 3D point. Also
              includes code for evaluating the corresponding jacobian
eucsbademo.h: function prototypes

calib.txt:    intrinsic calibration matrix for the employed camera

7cams.txt:    initial motion parameters for a test case involving 7 cameras
7pts.txt:     initial structure parameters for a test case involving 7 cameras

9cams.txt:    initial motion parameters for a test case involving 9 cameras
9pts.txt:     initial structure parameters for a test case involving 9 cameras

54cams.txt:    initial motion parameters for a test case involving 54 cameras
54pts.txt:     initial structure parameters for a test case involving 54 cameras

==================== COMPILING ====================
The demo program should be compiled during sba's compilation

==================== SAMPLE RUNS ====================
The command
  eucsbademo 7cams.txt 7pts.txt calib.txt
produces the following output:
  SBA using 465 3D pts, 7 frames and 1916 image projections

  Method BA_MOTSTRUCT, expert driver, analytic jacobian

  SBA returned 20 in 20 iter, reason 2, error 0.675396 [initial 19.0947]
  Elapsed time: 0.50 seconds, 500.00 msecs


For the 9 cameras case,
  eucsbademo 9cams.txt 9pts.txt calib.txt
produces
  SBA using 559 3D pts, 9 frames and 2422 image projections

  Method BA_MOTSTRUCT, expert driver, analytic jacobian

  SBA returned 22 in 22 iter, reason 2, error 0.619559 [initial 8.17604]
  Elapsed time: 0.69 seconds, 690.00 msecs


For the 54 cameras case,
  eucsbademo 54cams.txt 54pts.txt calib.txt
produces
  SBA using 5207 3D pts, 54 frames and 24609 image projections, 15999 variables

  Method BA_MOTSTRUCT, expert driver, analytic jacobian

  SBA returned 23 in 23 iter, reason 2, error 0.176473 [initial 2.14707]
  Elapsed time: 13.56 seconds, 13560.00 msecs


Send your comments/questions to lourakis@ics.forth.gr
