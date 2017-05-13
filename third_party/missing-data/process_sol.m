function [aff_err, shum_err, shum_aff_err, shum_stab] = process_sol(err, res, Merr, pts, INC, r)
% Take an initial solution and determine affine error, and errors after
% iterative refinement.
if err ~= -2
  aff_err = affine_error(res, pts, r);
  [shum_err, shum_res, shum_stab] = shum(Merr,INC,r,res,100);
  if shum_err ~= -2
    shum_aff_err = affine_error(shum_res, pts, r);
  else
    shum_aff_err = -2;
  end
else
  aff_err = -2;
  shum_err = -2;
  shum_aff_err = -2;
  shum_stab = 0;
end
