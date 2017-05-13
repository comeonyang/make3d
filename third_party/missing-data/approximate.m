function [error, Mapprox, stable1, stable2] = approximate(M, INC, r, NULLSPACE)

stable1 = 0;
stable2 = 0;
sumnonzeros = (sum(NULLSPACE' ~= 0)')';
samprows = find(0 ~= sumnonzeros);
nonsamprows = find(0 == sumnonzeros);

numcols = size(M,2);

[samp_approx,nonsampcols]=approx_matrix(M(samprows,:), INC(samprows,:), NULLSPACE(samprows,:),r,numcols);

% The rows of NULLSPACE now contain all the nullspaces of the crossproduct spaces.  
% Use the nullspace to approximate M, taking r least principal components.
if size(samp_approx,2) < r;
  error = -2;
% samp_approx = -2 means that for some reason approximation wasn't possible.
% or if < r columns, will not be able to extend the matrix correctly.
else
  sampcols = setdiff(1:numcols, nonsampcols);
  [Mapprox, stable1] ...
     = extend_matrix(M(:,sampcols), INC(:,sampcols), samp_approx, samprows, nonsamprows, r);
  [Mapprox_transpose, stable2] = extend_matrix(M', INC', Mapprox', sampcols, nonsampcols, r);
  Mapprox = Mapprox_transpose';

% These extensions are necessary because the nullspace might not allow us to compute
% every row of the approximating rank r linear space, and then this linear space might
% not allow us to fill in some columns, if they are missing too much data.

  error = sum(sum((M.*INC - Mapprox.*INC).^2));
end


