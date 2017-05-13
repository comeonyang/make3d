function S = random_scene(np)
% First we pick a scene with unit albedos, then scale by randomly chosen ones.
for i = 1:np
  US = [US, random_scene_normal];
end
albedos = rand(np,1);
S = [albedos, albedos, albedos]'.*US;
