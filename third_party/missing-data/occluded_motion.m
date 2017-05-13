function [M,INC, pts] = occluded_motion(fo,nframes,npoints,frot,ftrans)
% fo = fraction occluded.  We will make every point present for this fraction of
% the frames.
% nframes = number of frames.
% npoints = number of points.
% frot = fraction of 2pi degrees that rotation in depth occurs.
% ftrans = fraction of 1 for total in-plane translation
  
[M pts] = unoccluded_motion(nframes,npoints,frot,ftrans);
INC = motion_incidence(fo,nframes,npoints);
