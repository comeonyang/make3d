function [es, as, ss] = impute_compare2(frot,ftrans,sigma,num_frames,n)
% This will be a comparison with a different occlusion structure, of
% many points, clustered in groups of three or four, each appearing for very 
% few frames.

if num_frames == -1
  num_frames = 9;
end
ptframes = 5;
ptcols = 6;
numpts = ptcols*(1+num_frames-ptframes);

INC = [];
for i = 1:2:2*(num_frames+1-ptframes)
  INC = [INC, [zeros(i-1,ptcols); ones(2*ptframes,ptcols); ...
                zeros(2*num_frames-(i+2*ptframes-1),ptcols)]];
end

es = [];
as = [];
ss = [];
for i = 1:n
  [M, pts] = unoccluded_motion(num_frames, numpts, frot, ftrans);
%      pts = [rand(3,numpts);ones(1,numpts)];
%      M = rand(2*num_frames,4)*pts;
  ERR = randn(size(M)).*sigma;
  Merr = M + ERR;
  [e,a,s] = rankrsfm_aff_err(M, INC, ERR, pts);

Im = impute(Merr,INC);

  e2 = sum(sum((M.*INC - Im.*INC).^2));
  a2 = trans_affine_error(Im, pts);

numptsets = 1+num_frames-ptframes;
% This is the number of different rectangular blocks in INC.
midframe = 1+floor(numptsets/2);
% The frame where we start imputing.
halfrow = 1+2*(midframe-1);
halfcol = 1+ptcols*(midframe-1);
[lastrow,lastcol] = size(Merr);
Imhalf = impute(Merr(halfrow:lastrow, halfcol:lastcol), ...
                INC(halfrow:lastrow, halfcol:lastcol));
MerrImhalf = Merr;
MerrImhalf(halfrow:lastrow, halfcol:lastcol) = Imhalf;
INCImhalf = INC;
INCImhalf(halfrow:lastrow, halfcol:lastcol) = ...
  ones(1+lastrow-halfrow,1+lastcol-halfcol);
Im2 = impute_up(MerrImhalf,INCImhalf);

  e3 = sum(sum((M.*INC - Im2.*INC).^2));
  a3 = trans_affine_error(Im2, pts);

  es = [es;[e,e2,e3]];
  as = [as;[a,a2,a3]];
  ss = [ss;s];


end

global Mexternal INCexternal;
Mexternal = M;
INCexternal = INC;
