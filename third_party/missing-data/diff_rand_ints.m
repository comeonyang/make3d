function y = diff_rand_ints(ints, n,from,to)
%  Add n ints, randomly chosen between from and to inclusive to the list in ints.
%  The resulting vector of ints should have all different numbers.
if (to + 1 - from) < n + size(ints,2)
  error(sprintf('Not enough room for %d random numbers from %d to %d', n, from, to));
end
for i = 1:n
  x = random_int(from,to-size(ints,2));
  oldx = from - 1;
  while oldx < x
    temp = x;
    x = x + sum((oldx < ints) & (ints <= x));
    oldx = temp;
  end
  ints = [ints,x];
end
y = ints;