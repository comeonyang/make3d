function n = random_scene_normal()
p = [];
while p == []
  p = [(2*rand)-1, (2*rand)-1, rand]';
  if sum(p.*p) > 1 p = []; end
end
n=p./(sqrt(sum(p.*p)));