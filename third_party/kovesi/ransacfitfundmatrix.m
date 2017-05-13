% RANSACFITFUNDMATRIX - fits fundamental matrix using RANSAC
%
% Usage:   [F, inliers, NewDist, fail] = ransacfitfundmatrix(defaultPara, x1, x2, t, Depth1, Depth2, dist, feedback, disp, FlagDist)
%
% Arguments:
%          defaultPara - useful default parameters (like, camera intrinsic matrix)
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are normalised to that their
%                mean distance from the origin is sqrt(2).  The value of
%                t should be set relative to this, say in the range 
%                0.001 - 0.01  
%          Depth1/2 - depth imformation to support more accurate ransac
%          distrib  - initial distribution (default uniform dist)
%          feedback  - An optional flag 0/1. If set to one the trial count and the
%                 estimated total number of trials required is printed out at
%                 each step.  Defaults to 0.
%          disp  - if true, display the matches found when done.
%
%          FlagDist - if true, calculate the reprojection error
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%          NewDist - New Distribution after Ransac (Outliers have zero distribution)
%          fail    - true if Ransac fail to find any solution is not degenerated
%
% See Also: RANSAC, FUNDMATRIX

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

% February 2004  Original version
% August   2005  Distance error function changed to match changes in RANSAC
%
% additional parameter distrib is a vector of non-negative numbers that
% specifies a (not necessarily normalized) probability distribution over
% the different possible matches - passed to the ransac function to be
% used during the sampling procedure.
% (added by Jeff Michels)
% Additional NewDist estimated after Ransac is formed by using Depth
% information
% (added by Min Sun)

function [F, inliers, NewDist, fail] = ransacfitfundmatrix(defaultPara, x1, x2, t, Depth1, Depth2, distrib, feedback, disp, FlagDist, s)
% [F, inliers, fail] = ransacfitfundmatrix(x1, x2, t, feedback, distrib)

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    if nargin < 8
        feedback = 0;
        disp = 0;
        FlagDist = 0;
	    s = 8;  % Number of points needed to fit a fundamental matrix. Note that
		    % only 7 are needed but the function 'fundmatrix' only
		    % implements the 8-point solution.
    elseif nargin < 11
	    s = 8;  % Number of points needed to fit a fundamental matrix. Note that
		    % only 7 are needed but the function 'fundmatrix' only
		    % implements the 8-point solution.
    end
    
    [rows,npts] = size(x1);
    if rows~=2 & rows~=3
        error('x1 and x2 must have 2 or 3 rows');
    end
    
    if rows == 2    % Pad data with homogeneous scale factor of 1
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];   
    end
    
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.  Note that 'fundmatrix' will also call
    % 'normalise2dpts' but the code in 'ransac' that calls the distance
    % function will not - so it is best that we normalise beforehand.
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);

    
    fittingfn = @fundmatrix;
    distfn    = @funddist; % funciton handler below
    degenfn   = @isdegenerate; % function handler below
    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [F, inliers, NewDist, fail] = ...
        ransac(defaultPara, [x1; x2], [Depth1; Depth2], fittingfn, distfn, degenfn, s, t, distrib, [T1; T2], feedback, disp, FlagDist);

    % Now do a final least squares fit on the data points considered to
    % be inliers.
    F = fundmatrix(x1(:,inliers), x2(:,inliers));
    
    % Denormalise
    F = T2'*F*T1;
    
%--------------------------------------------------------------------------
% Function to evaluate the first order approximation of the geometric error
% (Sampson distance) of the fit of a fundamental matrix with respect to a
% set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
% 'Multiple View Geometry in Computer Vision', page 270.
%
% Note that this code allows for F being a cell array of fundamental matrices of
% which we have to pick the best one. (A 7 point solution can return up to 3
% solutions)
%
% Min add to calculate ReProjection Error from Depth information (4/22, 2007)

function [bestInliers, bestF, ReProjError] = funddist(F, x, t, defaultPara, Depth, T, FlagDist);
%[bestInliers, bestF] = funddist(F, x, t);

x1 = x(1:3,:);    % Extract x1 and x2 from x
x2 = x(4:6,:);
T1 = T(1:3,:);
T2 = T(4:6,:);

if iscell(F)  % We have several solutions each of which must be tested

  nF = length(F);   % Number of solutions to test
  bestF = F{1};     % Initial allocation of best solution
  ninliers = 0;     % Number of inliers

  for k = 1:nF
	t
    x2tFx1 = zeros(1,length(x1));
    for n = 1:length(x1)
      x2tFx1(n) = x2(:,n)'*F{k}*x1(:,n);
    end

    Fx1 = F{k}*x1;
    Ftx2 = F{k}'*x2;

    % Evaluate distances
    d =  x2tFx1.^2 ./ ...
      (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
	
    inliers = find(abs(d) < t);     % Indices of inlying points

    if length(inliers) > ninliers   % Record best solution
      ninliers = length(inliers);
      bestF = F{k};
      bestInliers = inliers;
    end
  end

else     % We just have one solution
  x2tFx1 = zeros(1,length(x1));
  for n = 1:length(x1)
    x2tFx1(n) = x2(:,n)'*F*x1(:,n);
  end

  Fx1 = F*x1;
  Ftx2 = F'*x2;

  % Evaluate distances
  d =  x2tFx1.^2 ./ ...
    (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
  figure(1);hist(d(d<t));
  bestInliers = find(abs(d) < t);     % Indices of inlying points
  bestF = F;                          % Copy F directly to bestF

end

% calculate ReProError ------------------------------------------------
if FlagDist
	Depth1 = Depth(1,:);
	Depth2 = Depth(2,:);
    X1 = inv(defaultPara.InrinsicK1)*inv(T1)*x1.*repmat(Depth1,3,1);
    X2 = inv(defaultPara.InrinsicK2)*inv(T2)*x2.*repmat(Depth2,3,1);
    E = (defaultPara.InrinsicK2'*T2'*bestF*T1*defaultPara.InrinsicK1);    
    [U S V] =svd(E);
    Rz_pos = [ [0 -1 0];...
               [1 0 0];...
               [0 0 1]];
    Rz_neg = Rz_pos';
    T_hat1 = U*Rz_pos*diag([1 1 0])*U';
    T1 = [-T_hat1(2,3);  T_hat1(1,3); -T_hat1(1,2)];
    R1 = U*Rz_pos'*V';
    T_hat2 = U*Rz_neg*diag([1 1 0])*U';
    T2 = [-T_hat2(2,3);  T_hat2(1,3); -T_hat2(1,2)];
    R2 = U*Rz_neg'*V';
    X1_2_1 = R1*X1;
    X1_2_2 = R2*X1;
    ops = sdpsettings('solver','sedumi','verbose',1);
    a = sdpvar(1,1);
    b = sdpvar(1,1);
    F = set(a>=0)+set(b>=0);
    sol = solvesdp(F,norm(X1_2_1(:)*a + b*repmat(T1, size(X1,2),1)- X2(:),1),ops);
    sol = solvesdp(F,norm(X1_2_2(:)*a + b*repmat(T2, size(X1,2),1)- X2(:),1),ops);
    a = double(a);
    b = double(b);
    ReProjError = sum(abs(X2_1*a - X1),1) + sum(abs(X1_2*b - X2),1);
    ReProjError(setdiff(1:size(ReProjError,2), bestInliers)) = Inf;
else
    ReProjError = [];
end 
% ---------------------------------------------------------------------`

%----------------------------------------------------------------------
% (Degenerate!) function to determine if a set of matched points will result
% in a degeneracy in the calculation of a fundamental matrix as needed by
% RANSAC.  This function assumes this cannot happen...
     
function r = isdegenerate(x)
    r = 0;    
    

