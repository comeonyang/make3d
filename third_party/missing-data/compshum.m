function [shum_err, shum_aff_err, shum_stab] = compshum(Merr, pts, INC, r, num_random_starts)
% Run Shum's alg. a number of times with different random starting points.
shum_stab = [];
shum_aff_err = [];
shum_err = [];

for i=1:num_random_starts
  [s_err, s_mat,s_stab] = shum(Merr,INC,r,0,100);
  if s_err ~= -2
    s_aff_err = affine_error(s_mat, pts, r);
  else
    s_aff_err = -2;
  end
  shum_stab = [ shum_stab, s_stab];
  shum_aff_err = [shum_aff_err,s_aff_err];
  shum_err = [shum_err, s_err];
end
