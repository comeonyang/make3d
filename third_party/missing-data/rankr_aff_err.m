function [err, aff_err, stable] = rankr_aff_err(M,INC,r,ERR,pts,n,iter)
% Run my method, and determine errors.

Merr = M + ERR;
INCsize = min(size(INC));
if INCsize < r
  err = -2;
  aff_err = -2;
else
  [rankr_err, rankr_res, stable] = rankr(Merr,INC,r,n, iter);
  if rankr_err ~= -2
    rankr_aff_err = affine_error(rankr_res, pts, r);
  else
    rankr_aff_err = -2;
  end
  err = rankr_err;
  aff_err = rankr_aff_err;
end
  
