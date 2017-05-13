function [M, INC,ERR] = shad_ints(fo, nframes, npoints,sigma)
% For now, I'll generate scenes with a surface normal drawn from 
% a uniform intensity facing the camera.  Lighting directions will
% be drawn from uniform density on some subset of the sphere to 
% generate the appropriate fo.

F = random_lights(nframes, fo);
P = random_scene(npoints);

ERR = randn(nframes,npoints).*sigma;
M = F*P;
INC = M+ERR > 0;
