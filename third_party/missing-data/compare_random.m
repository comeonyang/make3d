function y = compare_random(F,P,r,IP,sigma,it,n)

for j = 1:n
  U = rand(F,r);
  V = rand(P,r);
  M = U*V';
  ERR = randn(size(M)).*sigma;
  INC = rand(F,P)>IP;
  if rem(j,5)==0
    y = [y', compare(M,INC,r,ERR,it)']'
  else
    y = [y', compare(M,INC,r,ERR,it)']'
  end
end
