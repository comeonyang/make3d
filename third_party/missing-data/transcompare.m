function [err,afferr,stable] = transcompare(M1,M2,INC,r,ERR,pts,n,iter)
% run my algorithm (with an without iterative refinement) on two different 
% matrices.  M1 has no translation.  M2 has translation added, then approximately
% removed using the noisy, partial matrix.  So M1 is ground truth for both of them.
[e1,a1,s1] = rankr_aff_err(M1,INC,r,ERR,pts,n,iter);
[e2,a2,s2] = rankrsfm_aff_err(M2,INC,ERR,pts);
err = [e1, e2];
afferr = [a1,a2];
stable = [s1,s2];
