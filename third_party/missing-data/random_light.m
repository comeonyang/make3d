function l = random_light(fo)
if fo == .25
  dir = random_scene_normal;
elseif fo == .5
  dir = random_normal;
else
  error('Only implemented occlusion fractions of .25 and .5 for lighting.');
end
l = dir;
  
  