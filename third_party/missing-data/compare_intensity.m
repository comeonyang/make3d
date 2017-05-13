function y = compare_intensity(fo, nframes, npoints, sigma, s, n)
% fo = fraction of points we expect to be self-shadowed.  

num_random_starts = 5;

for j = 1:n
  [M,INC,ERR] = shad_ints(fo,nframes,npoints,sigma);
  if rem(j,5)==0
    y = [y', compare(M,INC,3,ERR,s,num_random_starts)']'
  else
    y = [y', compare(M,INC,3,ERR,s,num_random_starts)']';
  end
end
