function n = random_normal()
p = [];
while isempty(p)
  p = [(2*rand)-1, (2*rand)-1, (2*rand)-1]';
  if sum(p.*p) > 1 p = []; end
end
n=p./(sqrt(sum(p.*p)));