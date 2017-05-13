function [T, P] = tk_affine(M)
% This basically reimplements the Tomasi and Kanade algorithm for the case
% of affine motion, including allowing for translation.  We return the 
% motion matrix, with translation in the fourth column, and the pts, with 
% the 4th row = (1,1,1 ... 1,1).  The estimated matrix is T*P.

trans = -(sum(M'))/size(M,2);
transmat = repeat(trans',size(M,2));
M3 = M + transmat;

[u,s,v] = svd(M3);
vp = v';
u3 = u(:,1:3);
s3 = s(1:3,1:3);
s3sqrt = sqrt(s3);
vp3 = vp(1:3,:);

T = [u3*s3sqrt,-trans'];
P = [s3sqrt*vp3; ones(1,size(M,2))];

