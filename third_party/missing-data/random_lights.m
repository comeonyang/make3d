function L = random_lights(num, fo)
for i = 1:num
  L = [L, random_light(fo)];
end
L = L';