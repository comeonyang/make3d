function y = compare(M,INC,r,ERR,n,starts)
% run shum's method, and mine on the same matrix, and compare with ground truth.
% Gaussian error is added to the matrix, with sigma as given.
global Mexternal INCexternal

num_random_starts = starts;
num_data = num_random_starts+2+2;
% random starts + my method + my method & iterative, + ground truth + g.t. & it.

[M, INC, ERR] = rem_inc_els(M, INC, ERR, r);
Merr = M + ERR;
Mexternal = Merr;
INCexternal = INC;
INCsize = min(size(INC));
if INCsize < r
  y = repeat(-2,num_data);
else
  shum_err = [];
  for i=1:num_random_starts
    shum_err = [shum_err,shum(Merr,INC,r,0,100)];
  end
  [rankr_err, rankr_res] = rankr(Merr,INC,r,n);
  if rankr_err == -2
    rankr_shum_err = -2;
  else
    rankr_shum_err = shum(Merr,INC,r,rankr_res,100);
  end
%  if r==2
%    [rankall_err, rankall_res] = rank2(Merr,INC);
%  elseif r==3
%    [rankall_err, rankall_res] = rank3(Merr,INC);
%  else
%    rankall_err = -2;
%  end
%  if rankall_err == -2
%    rankall_shum_err = -2;
%  else
%    rankall_shum_err = shum(Merr,INC,r,rankall_res,100);
%  end
  ground_truth_err = sum(sum((M.*INC - Merr.*INC).^2));
  gt_shum_err = shum(Merr,INC,r,M,100);
  y=[shum_err,rankr_err,rankr_shum_err,ground_truth_err,gt_shum_err];
end