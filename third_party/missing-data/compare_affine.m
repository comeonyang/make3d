function [err, aff_err, stable] = compare_affine(M,INC,r,ERR,pts,n,starts, iter)
% run shum's method, and mine on the same matrix, and compare with ground truth.
% Gaussian error is added to the matrix, with sigma as given.
global Mexternal INCexternal

num_random_starts = starts;
num_data = num_random_starts+2+2;
% Each column of the results shows a different method; although the first
%   num_random_starts columns are all the same iterative method with different
%   random starting points.  Overall, they're organized as:
% random starts + my method + my method & iterative, + ground truth + g.t. & it.

[M, INC, ERR, pts] = rem_inc_els_pts(M, INC, ERR, pts, r);

occpercent = 1 - sum(sum(INC))/size(M,1)/size(M,2);

Merr = M + ERR;
Mexternal = Merr;
INCexternal = INC;
INCsize = min(size(INC));
if INCsize < r
  err = repeat(-2,num_data);
  aff_err = repeat(-2,num_data);
else
  shum_err = [];
  shum_aff_err = [];
  for i=1:num_random_starts
    [s_err, s_mat] = shum(Merr,INC,r,0,100);
    if s_err ~= -2
      s_aff_err = affine_error(s_mat, pts, r);
    else
      s_aff_err = -2;
    end
    shum_aff_err = [shum_aff_err,s_aff_err];
    shum_err = [shum_err, s_err];
  end
  [rankr_err, rankr_res, stable] = rankr(Merr,INC,r,n, iter);
  if rankr_err ~= -2
    rankr_aff_err = affine_error(rankr_res, pts, r);
  else
    rankr_aff_err = -2;
  end

  if rankr_err == -2
    rankr_shum_err = -2;
    rankr_shum_aff_err = -2;
  else
    [rankr_shum_err, rankr_shum_res] = shum(Merr,INC,r,rankr_res,100);
    if rankr_shum_err ~= -2
      rankr_shum_aff_err = affine_error(rankr_shum_res, pts, r);
    else
      rankr_shum_aff_err = -2;
    end
  end

  ground_truth_aff_err = affine_error(M,pts,r);
  ground_truth_err = sum(sum((M.*INC - Merr.*INC).^2));
  [gt_shum_err, gt_shum_res] = shum(Merr,INC,r,M,100);

  if gt_shum_err ~= -2
    gt_shum_aff_err = affine_error(gt_shum_res, pts, r);
  else
    gt_shum_aff_err = -2;
  end
  err=[shum_err,rankr_err,rankr_shum_err,ground_truth_err,gt_shum_err];
  aff_err=[shum_aff_err,rankr_aff_err,rankr_shum_aff_err,ground_truth_aff_err,gt_shum_aff_err];
end
