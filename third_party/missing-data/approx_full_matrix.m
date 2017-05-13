function [L, pts] = approx_full_matrix(M,r)
% approximate a matrix, M, that has no missing elements, with a rank r matrix
% this is essentially the Tomasi and Kanade algorithm.

[u,s,v] = svd(M);
vp = v';
% M = u*s*vp
ur = u(:,1:r);
sr = s(1:r,1:r);
vpr = vp(1:r,:);
L = ur*sr*vpr;
pts = [vpr; ones(1,size(vpr,2))];
