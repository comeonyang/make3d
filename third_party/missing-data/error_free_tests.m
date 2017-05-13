function y = error_free_tests(F,P,r,IP,n)
% Create a rank r FxP matrix by multiplying random Fxr and rxP matrices together.
% IP is fraction of data that will be missing from this matrix.  We check
% error produced by Shum algorithm n times, and return results.

% KNOWN PROBLEM: If incidence matrix is so sparse that there aren't enough
% values in each column and row something odd is happening.

global Mexternal INCexternal

for j = 1:n
  U = rand(F,r);
  V = rand(P,r);
  M = U*V';
  INC = rand(F,P)>IP;
  valid_cols = find(sum(INC)>r);
  valid_rows = find(sum(INC')>r);
  M = M(valid_rows,valid_cols);
Mexternal = M;
  INC = INC(valid_rows,valid_cols);
INCexternal = INC;
  error_shum = shum(M,INC,r)
  y_shum = [y_shum;error_shum];
  error_rank2 = rank2(M,INC)
  y_rank2 = [y_rank2;error_rank2];
end

y = [y_shum, y_rank2];