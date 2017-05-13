function l = occlusion_length(el,ml)
% generate a random length of occlusion.  el should be the expected
% length, ml should be the max length.  I'm going to make this very
% easy, by picking a uniform distribution of the max possible size.
if 2*el <= ml
  l = prob_round(2*el*rand);
else
  l = prob_round((2*el-ml)+2*(ml-el)*rand);
end
