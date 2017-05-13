function [M,INC, pts] = occluded_motion2(fo,nframes,npoints,frot,ftrans)
% fo = fraction occluded.  We will make every point present for this fraction of
% the frames.
% nframes = number of frames.
% npoints = number of points.
% frot = fraction of 2pi degrees that rotation in depth occurs.
% ftrans = fraction of 1 for total in-plane translation
  
% [M pts] = unoccluded_motion2(nframes,npoints,frot,ftrans);
[M pts] = arbitrary_motion(nframes,npoints,ftrans);  frot;
% ftrans isn't actually used here, 
INC = motion_incidence(fo,nframes,npoints);
