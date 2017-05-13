function [M, pts] = arbitrary_motion(nframes, npoints, trans)
% This will produce a completely arbitrary motion sequence, in which each frame has 
% nothing to do with the previous one.  

pts = rand(4,npoints);
pts(4,:) = ones(1,npoints);
% points have x,y,z coordinates, and a 1 in the 4th row for translation.

R = [];
for i = 1:nframes
  r1 = random_normal';
  % Produces a random unit vector.
  temp = random_normal';
  r2 = temp - sum(r1.*temp)*r1;
  R = [R;r1;r2];
end

if trans
  R = [R, 2*rand(2*nframes,1)-1];
else
  R = [R,zeros(2*nframes,1)];
end

M = R*pts;
