function [err, aff_err, stable] = rankrsfm_aff_err(M,INC,ERR,pts)
% Run my method, and determine errors.
Merr = M + ERR;
INCsize = min(size(INC));
if INCsize < 4
  err = -2;
  aff_err = -2;
  stable = 1;
else
%  [rankr_err, rankr_res, stable] = rankrsfm(Merr,INC);
  [rankr_err, rankr_res, stable] = rankrsfm_tpose(Merr,INC,1);
  if rankr_err ~= -2
    rankr_aff_err = trans_affine_error(rankr_res, pts);
  else
    rankr_aff_err = -2;
  end
  err = rankr_err;
  aff_err = rankr_aff_err;
end

% M - rankr_res
