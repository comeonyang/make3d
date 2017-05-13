function [M, pts] = unoccluded_motion(nframes,npoints,frot,ftrans)

% To generate this, we'll consider a very simple situation in which the motion
%   is constant, and the points are randomly chosen from within a cube.
%   the cube in between 0 and 1 in each dimension, so the image will be
%   roughly of size 1x1.
%   We assume that the key variable is the amount of rotation in depth, 
%   so we'll make the in-plane rotation constant so that it adds up to
%   ninety degrees, and the translation constant, 
%   in a random direction.

M = [];
translation_mag = ftrans/(nframes-1);
ip_rot = (pi/2)/(nframes-1);
op_rot = (2*pi*frot)/(nframes-1);

pts = rand(4,npoints);
pts(4,:) = ones(1,npoints);
% points have x,y,z coordinates, and a 1 in the 4th row for translation.


depth_axis = 2*pi*rand(1);
%depth_axis = 0;
translation_dir = 2*pi*rand(1);
translation_vec = translation_mag*[cos(translation_dir), sin(translation_dir),0];

rpts1 = rotate_z_axis(pts,depth_axis);


for fnum=0:nframes-1
  rpts2 = rotate_x_axis(rpts1,op_rot*fnum);
  rpts3 = rotate_z_axis(rpts2,-depth_axis);
  rpts4 = rotate_z_axis(rpts3,ip_rot*fnum);
  rpts5 = translate(rpts4,fnum*translation_vec);
  M = [M',rpts5(1:2,:)']';
end

