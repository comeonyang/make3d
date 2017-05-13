%% Copyright (c) 2001, NEC Research Institute Inc.  All rights reserved. 
%% 
%% Permission to use, copy, modify, and distribute this software and its
%% associated documentation for non-commercial purposes is hereby
%% granted, provided that the above copyright notice appears in all
%% copies, derivative works or modified versions of the software and any
%% portions thereof, and that both the copyright notice and this
%% permission notice appear in the documentation.  NEC Research Institute
%% Inc. shall be given a copy of any such derivative work or modified
%% version of the software and NEC Research Institute Inc. and its
%% affiliated companies (collectively referred to as NECI) shall be
%% granted permission to use, copy, modify and distribute the software
%% for internal use and research.  The name of NEC Research Institute
%% Inc. and its affiliated companies shall not be used in advertising or
%% publicity related to the distribution of the software, without the
%% prior written consent of NECI.  All copies, derivative works or
%% modified versions of the software shall be exported or reexported in
%% accordance with applicable laws and regulations relating to export
%% control.  This software is experimental.  NECI does not make any
%% representations regarding the suitability of this software for any
%% purpose and NECI will not support the software.  THE SOFTWARE IS
%% PROVIDED AS IS.  NECI DOES NOT MAKE ANY WARRANTIES EITHER EXPRESS OR
%% IMPLIED WITH REGARD TO THE SOFTWARE.  NECI ALSO DISCLAIMS ANY WARRANTY
%% THAT THE SOFTWARE IS FREE OF INFRINGEMENT OF ANY INTELLECTUAL PROPERTY
%% RIGHTS OF OTHERS.  NO OTHER LICENSE EXPRESS OR IMPLIED IS HEREBY
%% GRANTED. NECI SHALL NOT BE LIABLE FOR ANY DAMAGES, INCLUDING GENERAL,
%% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, ARISING OUT OF THE USE
%% OR INABILITY TO USE THE SOFTWARE.


% This README file will explain a little bit about my code for linear fitting
% with missing data.  This code addresses the problem of fitting a low rank
% approximation to a matrix in which some of the elements are unknown.  That is,
% this algorithm attacks the problem that is easily solved by SVD when there
% is no missing data.  This directory also contains code that is specialized to 
% solving structure from motion (SFM) problems.  A full description of the algorithms 
% can be found in ``Linear Fitting with Missing Data for Structure-From-Motion,'' 
% Computer Vision and Image Understanding, 82:57-81, (2001), D. Jacobs.  
% An earlier version of this work, with earlier code, was released in '97 and
% presented in CVPR '97.  That code did not have the specializations for SFM,
% and did not work on motion sequences with translation (which requires a bit
% more than a generic low rank approximation).

% This directory contains not just the code for my algorithm, but also the 
% code used to generate all the tests described in the CVIU paper, and
% some additional tests as well.  In particular, there is also code that 
% implements an iterative method due to Wiberg, which is described in a paper
% by Shum, Ikeuchi and Reddy in PAMI, Sept. 1995.  And there is code to 
% generate synthetic motion sequences, occlusions and intensity images.  I am
% only trying to introduce users to the central code in this file; if you want
% to try to use the code developed for my experiments you'll have to nose around
% these files.

% Note that my earlier code (rankr.m) does have a lot of little hacks that are not
% well documented.  For the most part these concern choosing column triples
% to estimate the nullspace, and methods for extending an initial stable solution 
% when the stable solution only covers part of the matrix.  This typically occurs
% when the amount of missing data is large.

% COMPATIBILITY NOTE: This code has been upgraded to run in MATLAB V.  As
% a consequence, it won't run in MATLAB IV.  To get it to run in MATLAB IV
% you must:

% -- Remove the call to "logical" in whole_invert.  

% On the other hand, I have not completely tested the port to MATLAB V.
% However, the main functions, shown below, have been tested as shown.

% This file is also a Matlab M-file, and can be run to see a quick demo of
% some of the main functions.

% Although I don't promise to support this code, or upgrade, any questions
% or comments can reach me at: djacobs@cs.umd.edu.

% First, I'll demonstrate the generic code on a special case for which it is 
% particularly good.  

M = 100.*(rand(6,3)*(rand(3,9)))
INC = [1,1,1,1,1,1,0,0,0; ...
       1,1,1,1,1,1,0,0,0; ...
       1,1,1,0,0,0,1,1,1; ...
       1,1,1,0,0,0,1,1,1; ...
       0,0,0,1,1,1,1,1,1; ...
       0,0,0,1,1,1,1,1,1]

% We've just generated a small, 6x9 matrix that is rank 3.  
% INC provides a pattern of missing data for this matrix, with 
% 1 indicating the data is present, and 0 that it is missing.
% None of the algorithms will look at the data associated with
% 0s, except the code to evaluate the results.  This matrix has
% the interesting property that there is no 4x4 subblock that is missing 
% only one element, so the Tomasi and Kanade approach of "hallucinating"
% the values of missing elements will not work.  To apply my algorithm
% to this data, one can type:

disp('Our solution is:');
[e,a] = rankr(M,INC,3,1000)

% This finds a rank 3 approximation to M, with the values in INC missing.
% The fourth parameter, 1000, says that 1000 triples of columns are the maximum
% number that will be used in estimating the nullspace; if we don't have enough 
% information after that, we give up (see the paper if this doesn't make
% any sense).  For this particular example, it would make more sense to just try 
% all triples of columns, instead of a random sample, but for big matrices that 
% won't work.  However, because of this randomness, there's a small chance that the
% algorithm won't find the right answer.
% The value e that is returned is the sum of square differences between the
% elements that are present in the matrix, and the corresponding values in the
% approximating matrix.  a is the approximating matrix.

% This matrix also has too many missing elements for the iterative method of
% Shum et al to work.  Here's how we call it:

disp('The iterative method gives:');
[e,a] = shum(M,INC,3,0,100)

% Here, e = -2 signals that the algorithm couldn't find a solution.  
% The parameter 3 means find a rank 3 approximation, the parameter 0
% means pick a random matrix as the starting point (we could use a starting guess
% here), and 100 means that if the algorithm doesn't converge after 100 iterations,
% stop anyway.

S = random_motion(4);
P = [rand(3,12); ones(1,12)];
% S contains random motions for 4 frames.  P contains 12 random points.
% The ones are there so S*P gives the image points.
disp('M are image points from a random motion');
M = S*P
INC = [1,1,1,1,1,1,1,1,1,0,0,0;...
      1,1,1,1,1,1,1,1,1,0,0,0;...
      1,1,1,1,1,1,0,0,0,1,1,1;...
      1,1,1,1,1,1,0,0,0,1,1,1;...
      1,1,1,0,0,0,1,1,1,1,1,1;...
      1,1,1,0,0,0,1,1,1,1,1,1;...
      0,0,0,1,1,1,1,1,1,1,1,1;...
      0,0,0,1,1,1,1,1,1,1,1,1]
% This also represents a pattern of missing data that our algorithm can handle that 
% is tricky for other approaches.  See the paper for details.

[e,a,s] = rankrsfm(M,INC)


% Now, let's try an example that all methods work on.
M = rand(20,3)*rand(3,20);
disp('We set M to be a rank 3, 20x20 matrix.');
INC = rand(20,20)>.25;
disp('We set INC to be a random matrix, with each element 1 with prob. .75.');

disp('Our solution has error of:');
[e,a] = rankr(M,INC,3,1000);
disp(e);

disp('If we use it as the starting point for an iterative method, we get error of:');
[e1,a1] = shum(M,INC,3,a,100);
disp(e1);

disp('If we run the iterative method with a random starting point, we get error of:');
[e2,a2] = shum(M,INC,3,0,100);
disp(e2);

% Although results are randomized, probably all three approaches gave a good result
% here.

% There is also a lot of code to run experiments, generate synthetic test data, and
% compare methods.  If you look at compare_motion, you'll get a starting point to
% the code I used to test the generic method for the CVPR '97 paper.  compare_transmot
% contains code used to run experiments on the new SFM code.

% Good luck
% -- David Jacobs
