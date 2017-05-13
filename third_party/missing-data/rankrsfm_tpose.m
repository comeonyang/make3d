function [error, Mapprox, stable] = rankrsfm_tpose(M,INC,trans)
% Compare to rankrsfm.  This will work on the transpose of the matrix.
% trans is a switch that allows the same routine to either allow for 
% translation or not.

if ~ trans
  r = 3;
else 
  r = 4;
end
M = M';
INC = INC';
MAXNULLSIZE = 10;
[numrows,numcols] = size(M);
NULLSPACE = [];
frame_pairs = random_frame_pairs(numcols/2);
% gives a 2 x (numframes choose 2) matrix with all pairs of frames, in random order.

for pair = frame_pairs
  if size(NULLSPACE,2) >= numcols*MAXNULLSIZE, break, end

  [submatrix,subrows] = fp_submatrix(pair,M',INC');
  submatrix = submatrix';
  % I have to transpose everything, since this routine is written to work with
  % untransposed matrix.
  num_subpoints = size(submatrix,1);
  if num_subpoints > r
    if trans
      submatrix = [ones(num_subpoints,1), submatrix];
    end
    subnull = nulleps(submatrix,.001,r);
    if size(subnull,2) == num_subpoints - r
       % Since submatrix has 5+ rows and 5 columns, if it has rank 4, subnull is 
       % number rows - 4.
       % If subnull is bigger, submatrix is rank deficient, and results unstable.
       nullTEMP = zeros(numrows, size(subnull,2));
       nullTEMP(subrows,:) = subnull;
       NULLSPACE = [NULLSPACE, nullTEMP];
    end
  end
end

if size(NULLSPACE) == [0 0]
  error = -2;
  Mapprox = [];
else
  [error, Mapprox, stable1, stable2] = approximate(M,INC,r,NULLSPACE);
  % Now we should be sure that the nullspace of NULLSPACE contains a column with all ones.
  if error ~= -2
    error = sum(sum((M.*INC - Mapprox.*INC).^2));
    Mapprox = Mapprox';
  end
  stable = stable1 & stable2;
end

if error == -2
  stable = 0;
end




  

