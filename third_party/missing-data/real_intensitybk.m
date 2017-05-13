function [err, aff_err] = real_intensity(s,n)
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

%  scaled_corner_data;
%  M = M(:,[1:100,201:400]);

  test2_data;
  M = test2_im_data;

  [nr,m] = size(M);
  ERR = zeros(nr,m);

%  for i = 1:100
%    norms(:,i) = [1 0 0 1]';
%  end
%  for k = 1:10
%    for i = 1:10
%      norms(:,100 + (k-1)*20 + i) = [0 1 0 1]';
%    end
%    for i = 11:20
%      norms(:,100 + (k-1)*20 + i) = [0 0 1 1]';
%    end
%  end

  for i = 1:432
    norms(:,i) = [1 0 0 1]';
  end

for j = 1:n

%  INC = rand(nr,m) > fo; 
  INC = (test2_im_data > -.8);

  if rem(j,1000)==0
    [e, a] = compare_affine(M,INC,3,ERR,norms,s,num_random_starts);
    err = [err', e']'
    aff_err = [aff_err', a']'
  else
    [e, a] = compare_affine(M,INC,3,ERR,norms,s,num_random_starts);
    err = [err', e']';
    aff_err = [aff_err', a']';
  end
end

