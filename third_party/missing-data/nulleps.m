function N = nulleps(M,tol,maxrank)
% Find the nullspace of M.  This is the regular nullspace, augmented by 
% vectors that aren't really in the nullspace, but have low singular
% values associated with them.  tol is the threshold on these singular values.

% maxrank indicates the maximum rank that M can be expected to have.  All but the 
% first MAXRANK singular values can be assumed to be due to noise.  Ideally, this 
% should be treated in a unified way with the tolerance, and one should detect 
% whether results seem unstable when the MAXRANK singular value is near to the 
% MAXRANK+1.
if maxrank == -1
  maxrank = min(size(M));
end

[u,s,v] = svd(M);
sigsvs = min(maxrank, sum(diag(s)>tol));
N=u(:,sigsvs+1:size(u,2));
