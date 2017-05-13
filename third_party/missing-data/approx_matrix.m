function [Mapprox,misscols] = approx_matrix(M,INC,NULLSPACE,r,numcols)
Mapprox = [];
misscols = [];

%[A,S,U] = svd(NULLSPACE',0);
[U,S,V] = svd(NULLSPACE);

Unumcols = size(U,2);
%SINGULARVALUES = diag(S)'
if svd_suff_data(S,r)
  RSPACE = U(:,Unumcols+1-r:Unumcols);
% RSPACE has r columns, which span the r-D space that gives the best linear
% surface to approximate M.
  for i = 1:numcols
    INCcol = INC(:,i);
    Mcol = M(:,i).*INCcol;
    INCmat = repeat(INCcol, r);
    RSPACEcols = RSPACE.*INCmat;
    % this is just RSPACE with all the rows missing that are also missing in Mcol.
    if rank(RSPACEcols) == r
      Mapprox = [Mapprox, RSPACE*(RSPACEcols\Mcol)];
    else
      misscols = [misscols, i];
    end
  end
else
  Mapprox = -2;
end
if isempty(Mapprox) Mapprox = -2; end