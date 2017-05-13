% Demonstration of feature matching via simple correlation, and then using
% RANSAC to estimate the fundamental matrix and at the same time identify
% (mostly) inlying matches
%
% Usage:  testfund              - Demonstrates fundamental matrix calculation
%                                 on two default images
%         testfund(im1,im2)     - Computes fundamental matrix on two supplied images
%
% Edit code as necessary to tweak parameters

% Copyright (c) 2004-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% February 2004
% August   2005 Octave compatibility

function [Ftrans,correlations]=testfund(im1,im2)
    
    if nargin == 0
	im1 = imread('im02.jpg');
	im2 = imread('im03.jpg');
    end
    correlations = [];
    v = version; Octave=0;% Crude Octave test        
    thresh = 500;   % Harris corner threshold
    nonmaxrad = 3;  % Non-maximal suppression radius
    dmax = 100;      % Maximum search distance for matching
    w = 11;         % Window size for correlation matching
    
    % Find Harris corners in image1 and image2
    [cim1, r1, c1] = harris(im1, 1, thresh, 3); clear cim1;
    show(im1,1), hold on, plot(c1,r1,'r+');

    [cim2, r2, c2] = harris(im2, 1, thresh, 3); clear cim2;
    show(im2,2), hold on, plot(c2,r2,'r+');
    drawnow

    correlation = 0;  % Change this between 1 or 0 to switch between the two
                      % matching functions below
    
    if correlation  % Use normalised correlation matching
	[m1,m2] = matchbycorrelation(im1, [r1';c1'], im2, [r2';c2'], w, dmax);
	
    else            % Use monogenic phase matching
	nscale = 1;
	minWaveLength = 10;
	mult = 4;
	sigmaOnf = .2;
	[m1,m2] = matchbymonogenicphase(im1, [r1';c1'], im2, [r2';c2'], w, dmax,...
					nscale, minWaveLength, mult, sigmaOnf);
    end    
    
    % Display putative matches
    show(im1,3), set(3,'name','Putative matches')
    if Octave, figure(1); title('Putative matches'), axis('equal'), end        
    for n = 1:length(m1);
	line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)])
    end

    % Assemble homogeneous feature coordinates for fitting of the
    % fundamental matrix, note that [x,y] corresponds to [col, row]
    x1 = [m1(2,:); m1(1,:); ones(1,length(m1))];
    x2 = [m2(2,:); m2(1,:); ones(1,length(m1))];    
    
    t = .002;  % Distance threshold for deciding outliers
    
    % Change the commenting on the lines below to switch between the use
    % of 7 or 8 point fundamental matrix solutions, or affine fundamental
    % matrix solution.
%   [F, inliers] = ransacfitfundmatrix7(x1, x2, t, 1);    
   [F, inliers] = ransacfitfundmatrix(x1, x2, t, 1);
%   [F, inliers] = ransacfitaffinefund(x1, x2, t, 1);    
    fprintf('Number of inliers was %d (%d%%) \n', ...
	    length(inliers),round(100*length(inliers)/length(m1)))

    [Ftrans, transinliers] = ransacfittransfundmatrix(x1, x2, t, 1);
    correlations = [correlations;[m1(:,transinliers)',m2(:,transinliers)']];
    [U,S,V]=svd(Ftrans);
    epipole=[V(2,3)/V(3,3);V(1,3)/V(3,3)]
    fprintf('Number of inliers was %d (%d%%) \n', ...
	    length(inliers),round(100*length(transinliers)/length(m1)))
    fprintf('Number of putative matches was %d \n', length(m1))        
    
    % Display both images overlayed with inlying matched feature points
    
    if Octave
      figure(4); title('Inlying matches'), axis('equal'),
    else
      show(im1,4), set(4,'name','Inlying matches'), hold on
    end
    plot(m1(2,inliers),m1(1,inliers),'c+');
    %    plot(m2(2,inliers),m2(1,inliers),'g+');

    for n = inliers
      line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end

    show(im2,5), set(5,'name','Inlying matches'), hold on
    %    plot(m1(2,inliers),m1(1,inliers),'c+');
    plot(m2(2,inliers),m2(1,inliers),'g+');

    for n = inliers
      line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end


    if Octave
      figure(4); title('Inlying matches'), axis('equal'),
    else
      show(im1,6), set(6,'name','Translational Inlying matches'), hold on
    end
    plot(m1(2,transinliers),m1(1,transinliers),'c+');
    %    plot(m2(2,inliers),m2(1,inliers),'g+');

    for n = transinliers
      line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end

    show(im2,7), set(7,'name','translational Inlying matches'), hold on
    %    plot(m1(2,inliers),m1(1,inliers),'c+');
    plot(m2(2,transinliers),m2(1,transinliers),'g+');

    for n = transinliers
      line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end
  
  % determine which picture is closer to the epipole
  p1 = [r1(transinliers)-epipole(1),c1(transinliers)-epipole(2)]';
  p2 = [r2(transinliers)-epipole(1),c2(transinliers)-epipole(2)]';
  dist1 = sum((sum(p1.^2)).^.5);
  dist2 = sum((sum(p2.^2)).^.5);
  cImage = (dist1>dist2)+1;
  
  DualMatches=guide([r1,c1]', [r2,c2]', epipole, 17, .85, .2,150,cImage,im1,im2 );
%  correlations = [correlations;[r1(DualMatches(:,1)),c1(DualMatches(:,1)),r2(DualMatches(:,2)),c2(DualMatches(:,2))]];
  correlations = [r1(DualMatches(:,1)),c1(DualMatches(:,1)),r2(DualMatches(:,2)),c2(DualMatches(:,2))];
  fprintf('Number of guided matches was %d \n', size(DualMatches,1))        
    
  show(im1,9), set(9,'name','Dual matches'), hold on
      plot(c1(DualMatches(:,1)),r1(DualMatches(:,1)),'g+');
  %    plot(c2(DualMatches(:,2)),r2(DualMatches(:,2)),'g+');    

  %plot(epipole(2),epipole(1),'r*');    
  for n = 1:1:size(DualMatches,1)
    line([c1(DualMatches(n,1)) c2(DualMatches(n,2))], [r1(DualMatches(n,1)) r2(DualMatches(n,2))],'color',[0 0 1])
  end

  show(im2,8), set(8,'name','Dual matches'), hold on
  %    plot(c1(DualMatches(:,1)),r1(DualMatches(:,1)),'c+');
      plot(c2(DualMatches(:,2)),r2(DualMatches(:,2)),'g+');    

  %plot(epipole(2),epipole(1),'r*');    
  for n = 1:1:size(DualMatches,1)
    line([c1(DualMatches(n,1)) c2(DualMatches(n,2))], [r1(DualMatches(n,1)) r2(DualMatches(n,2))],'color',[0 0 1])
  end

  show(im1,11), set(11,'name','All matches'), hold on
      plot(correlations(:,2),correlations(:,1),'g+');
  %    plot(c2(DualMatches(:,2)),r2(DualMatches(:,2)),'g+');    

  %plot(epipole(2),epipole(1),'r*');    
  for n = 1:1:size(correlations,1)
    line([correlations(n,2) correlations(n,4)], [correlations(n,1) correlations(n,3)],'color',[0 0 1])
  end

  show(im2,10), set(10,'name','All matches'), hold on
  %    plot(c1(DualMatches(:,1)),r1(DualMatches(:,1)),'c+');
      plot(correlations(:,4),correlations(:,3),'g+');    

  %plot(epipole(2),epipole(1),'r*');    
  for n = 1:1:size(correlations,1)
    line([correlations(n,2) correlations(n,4)], [correlations(n,1) correlations(n,3)],'color',[0 0 1])
  end
  n