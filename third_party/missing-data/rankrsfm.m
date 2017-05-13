function [error, Mapprox, stable] = rankrsfm(M,INC)
% rankr is a generic routine for fitting linear surfaces.  This routine is 
% designed the case of structure from motion, with translation.  This differs
% from the generic case in two ways.  First, translation implies that M has 
% rank r+1 (typically 4, unless motion is somehow limited; in this routine
% we will assume r=3, since this influences column selection strategies), and that 
% one of the basis vectors for this space is a row vector of (1 1 1 ... 1).
% Second, in structure from motion there must be a specific pattern of 
% missing data, since the x and y coordinates of a point are either both present, 
% or both missing.

% For the second reason, we will consider every pair of frames, rather than 
% looking for larger combinations of points.  

% We then have two possible strategies, depending on whether we work with the
% matrix or its transpose.

MAXNULLSIZE = 10;
[numrows,numcols] = size(M);
NULLSPACE = [];
frame_pairs = random_frame_pairs(numrows/2);
% gives a 2 x (numframes choose 2) matrix with all pairs of frames, in random order.

for pair = frame_pairs
  if size(NULLSPACE,2) >= numrows*MAXNULLSIZE, break, end

  [submatrix,subcols] = fp_submatrix(pair,M,INC);
  subrows = [2*pair(1,1)-1, 2*pair(1,1), 2*pair(2,1)-1, 2*pair(2,1)];
  % Right now, I'm only doing nullspace of matrix, not it's transpose.
  num_subpoints = size(submatrix,2);
  if num_subpoints > 3
    submatrix = [ones(1,num_subpoints); submatrix];
    subnull = nulleps(submatrix,.001,4);
    if size(subnull,2) == 1
       % Since submatrix has 5 rows and 4+ columns, if it has rank 4, subnull is 1 col.
       % If subnull is bigger, submatrix is rank deficient, and results unstable.
       nullTEMP = zeros(numrows+1, 1);
       nullTEMP([1,1+subrows],:) = subnull;
       NULLSPACE = [NULLSPACE, nullTEMP];
    end
  end
end

%NULLSPACE
%error = NULLSPACE;
%return;

% size(NULLSPACE)
% rank(NULLSPACE)
% global NULLSPACEexternal
% NULLSPACEexternal = NULLSPACE;

if size(NULLSPACE) == [0 0]
  error = -2;
  Mapprox = [];
else
  [error, Mapprox, stable1, stable2] ...
    = approximate([ones(1,numcols);M],[ones(1,numcols);INC],4,NULLSPACE);
  % This works, but it's not clear that the best use of knowledge that there's
  % translation is to just add a row of 1's.  One possibility is to make
  % this row have large magnitude, so that it will be very well approximated
  % (since it has no noise).  Another is to figure out a correct method.
  if error ~= -2
    Mapprox = Mapprox(2:numrows+1,:);
    error = sum(sum((M.*INC - Mapprox.*INC).^2));
  end
  % In determining error, error in this row shouldn't be counted.
  stable = stable1 & stable2
end

if error == -2
  stable = 0;
end
