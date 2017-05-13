function [toterr,totafferr,totstable] = comp_transmot(fo,nframes,npoints,frot,ftrans,sigma,n, real)
% I will compare my methods to see how they work when there is translation.
%   The parameters to this function are just like compare_motions.  I also allow 
%   a translation to be added to the motions (I don't think this makes any difference.
%   what's important is whether translation is known or not.  Whether there actually
%   is any translation or not, and what it's magnitude is, shouldn't matter).
%   ftrans indicates magnitude of translation (direction will be constant).  Reasonable 
%     ftrans is probably about .5.
%
% I will compare two possibilities.  One, in which there is no translation
%   (this is like the old experiments).  
%   Two, in which translation is added in to a motion that is otherwise the same, 
%   and the new method is used to estimate the rank 4 matrix with one
%   vector constrained to be (1 1 1  ... 1).
%
% In the first case, I will run SHUM with five random starting points; I will run 
% my new method with and without iterative refinement, my old method with and 
% without refinement (if this is clearly inferior, I may eliminate it) 
% and I'll take ground truth with and without iterative refinement.  
% In the second case, I'll just run the new method.

r = 3;
num_data = 12;
coltrips = 1000;

for j = 1:n
  if ~ isempty(real)
    INC = motion_incidence(fo,nframes,npoints);
    Dtj = remove_translations(real, INC);
    [M, pts] = approx_full_matrix(Dtj,3);
    ERR = Dtj - M;
    Mt = real - ERR;
  end
  if isempty(real) [M,INC,pts] = occluded_motion(fo,nframes,npoints,frot,0); end
  if isempty(real) ERR = randn(size(M)).*sigma; end
  [M, INC, ERR, pts] = rem_inc_els_pts(M, INC, ERR, pts, 4);
%  occpercent = 1 - sum(sum(INC))/size(M,1)/size(M,2);
  INCsize = min(size(INC));
  if INCsize < 4
    err = repeat(-1,num_data);
    aff_err = repeat(-1,num_data);
    stable = repeat(-1,num_data);
  else
    Mt = add_translations(M,ftrans);
    Merr = M + ERR;
    [shum_err, shum_aff_err, shum_stab] = compshum(Merr, pts, INC, r, 5);
    [old_err, old_res, old_stable] = rankr(Merr,INC,3,coltrips,j);
    [old_aff_err, old_shum_err, old_shum_aff_err, old_shum_stab] ...
       = process_sol(old_err, old_res, Merr, pts, INC, r);
    [rankr_err, rankr_res, rankr_stable] = rankrsfm_tpose(Merr,INC,0);
    % 0 means no translations.
    [rankr_aff_err, rankr_shum_err, rankr_shum_aff_err, rankr_shum_stab] ...
       = process_sol(rankr_err, rankr_res, Merr, pts, INC, r);
    [rankr_trans_err, rankr_trans_aff_err, trans_stable] = rankrsfm_aff_err(Mt,INC,ERR,pts);
    ground_truth_aff_err = affine_error(M,pts,r);
    ground_truth_err = sum(sum((M.*INC - Merr.*INC).^2));
    [gt_shum_err, gt_shum_res, gt_shum_stab] = shum(Merr,INC,r,M,100);
    if gt_shum_err ~= -2
      gt_shum_aff_err = affine_error(gt_shum_res, pts, r);
    else
      gt_shum_aff_err = -2;
    end

    err=[shum_err,old_err,old_shum_err,rankr_err,rankr_shum_err,rankr_trans_err, ...
         ground_truth_err,gt_shum_err];
    afferr=[shum_aff_err,old_aff_err,old_shum_aff_err,rankr_aff_err,rankr_shum_aff_err,...
               rankr_trans_aff_err,ground_truth_aff_err,gt_shum_aff_err];
    stable = [shum_stab, old_stable, old_shum_stab, rankr_stable, rankr_shum_stab, ...
               trans_stable, gt_shum_stab];
  end
  toterr = [toterr', err']';
  totafferr = [totafferr', afferr']';
  totstable = [totstable', stable']';
  if j < 11 | j == 20 | j == 30 | rem(j,100) == 0
    toterr
    totafferr
    totstable
  end
end
