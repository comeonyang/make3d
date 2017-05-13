function [M, pts] = unoccluded2_motion(nframes,npoints,frot,ftrans)
% This is an alternate motion generator.  The previous one had the property that the 
% if, for each frame, one takes the third row of the rotation matrix, these all lie in 
% a plane.

% To avoid this, we'll add a little torque to the old method.  We'll rotate out of plane by 
% the same magnitude in another direction, but this time with angular acceleration.
M = [];
translation_mag = ftrans/(nframes-1);
ip_rot = (pi/2)/(nframes-1);
op_rot = (2*pi*frot)/(nframes-1);

tor_rots = (0:nframes-1).*(1:nframes)./2;
tor_rots = (2*pi*frot).*tor_rots./sum(tor_rots);
tor_axis = 2*pi*rand(1);

pts = rand(4,npoints);
pts(4,:) = ones(1,npoints);
% points have x,y,z coordinates, and a 1 in the 4th row for translation.

depth_axis = 2*pi*rand(1);
translation_dir = 2*pi*rand(1);
translation_vec = translation_mag*[cos(translation_dir), sin(translation_dir),0];

rpts1 = rotate_z_axis(pts,depth_axis);

for fnum=0:nframes-1
  tor_rot = tor_rots(fnum+1);
  rpts2 = rotate_x_axis(rpts1,op_rot*fnum);
  rpts3 = rotate_z_axis(rpts2,-depth_axis);

  rpts3a = rotate_z_axis(rpts3,tor_axis);
  rpts3b = rotate_x_axis(rpts3a,tor_rot);
  rpts3c = rotate_z_axis(rpts3b,-tor_axis);

  rpts4 = rotate_z_axis(rpts3c,ip_rot*fnum);
  rpts5 = translate(rpts4,fnum*translation_vec);
  M = [M',rpts5(1:2,:)']';
end


