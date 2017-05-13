function [err, aff_err, stab] = compare_motion(fo,nframes,npoints,frot,sigma,s,n)
% fo = expected fraction of frames points are occluded.  The actual number
%   of frames a point appear in is generated randomly for each point.
% nframes = number frames
% npoints = number points
% frot = fraction of 2*pi for rotation in depth.
% sigma = magnitude of Gaussian error.  Pt. coordinates typically range
%   between -1 and 1, so a reasonable setting might be .01, or .005.
% s = number of samples to draw in my method.
% n = number times to repeat exp.

% stab tells, for each iteration, whether the rankr result was stable.

% there are a number of other parameters set in these experiments.
% In motion generation, there is amount of translation (0) and in-plane rotation.
%   in shum there are the number of iterations to try.

stab = [];
err = [];
aff_err = [];
num_random_starts = 5;

for j = 1:n
  [M,INC,pts] = occluded_motion(fo,nframes,npoints,frot,0);
  ERR = randn(size(M)).*sigma;
  if rem(j,1000)==0
    [e, a, stable] = compare_affine(M,INC,3,ERR,pts,s,num_random_starts, j);
    stab = [stab', stable']'
    err = [err', e']'
    aff_err = [aff_err', a']'
  else
    [e, a, stable] = compare_affine(M,INC,3,ERR,pts,s,num_random_starts, j);
    stab = [stab', stable']';
    err = [err', e']';
    aff_err = [aff_err', a']';
  end
end

