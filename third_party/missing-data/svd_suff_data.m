function y = svd_suff_data(S,r)
% S is the singular value part of the svd of the nullspaces of the column
% r-tuples.  We'll want to be able to take the r least significant columns
% of U.  This is right because the columns of U should span the whole space
% that M's columns might span.  That is, M is FxP.  The columns of U should
% span the F-dimensional Euclidean space, since U is FxF.  However, we want 
% to make sure that the F-r-1'th singular value of S isn't tiny.  If it is,
% our answer is totally unreliable, because the nullspaces of the column 
% r-tuples don't have sufficient rank.  If this happens, it means that the 
% intersection of the column cross-product spaces is bigger than r-dimensional,
% and randomly choosing an r-dimensional subspace of that isn't likely to
% give the right answer.
Snumcols = size(S,2);
Snumrows = size(S,1);
if (Snumrows == 0 | Snumcols + r < Snumrows)
  y = 0;
else
  y = S(Snumrows-r,Snumrows-r)>.001;
end
