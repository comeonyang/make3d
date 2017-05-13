function [err, aff_err, stable] = real_intensity(s)
% fo = fraction of frames points are occluded.  
% nframes = number frames
% npoints = number points
% frot = fraction of 2*pi for rotation in depth.
% sigma = magnitude of Gaussian error.  Pt. coordinates typically range
%   between -1 and 1, so a reasonable setting might be .01, or .005.
% s = number of samples to draw in my method.
% n = number times to repeat exp.

% there are a number of other parameters set in these experiments.
% In motion generation, there is amount of translation (0) and in-plane rotation.
%   in shum there are the number of iterations to try.

num_random_starts = 5;

test2_data;
M = test2_im_data;

ERR = zeros(size(M));


% norms is meaningless here, since we do not have ground truth
% geometric data (normals) for this image.  just make norms something,
%so the affine_error computation does not generate an error.

for i = 1:432
  norms(:,i) = [1 0 0 1]';
end

INC = (test2_im_data > .2);

[err, aff_err, stable] = compare_affine(M,INC,3,ERR,norms,s,num_random_starts,1);


