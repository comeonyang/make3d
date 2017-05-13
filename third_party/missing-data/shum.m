function [error, Wapp,stable] = shum(W,INC,r,INIT,num_it)
% shum produces an approximation to the W matrix, using method of Shum, Ikeuchi and Reddy
% W is a matrix of data.  INC is an incidence matrix of the same dimension.  If
%   INC is 1, the corresponding value in W is a real value.  Otherwise, it is really
%   missing data, and that entry in W should be ignored.
% r is the rank of the matrix we want to approximate W with.  
% INIT is an initial guess about the matrix, of the same size as W.  If INIT is 0,
%   we initialize with a random matrix, with elements from 0 to 1.
% We repeat for either num_it iterations, or till method converges to nearly
%   0 error.  
% We return the matrix, and the residual error.  By convention, if the method
% doesn't converge to 10 decimal places we return the negation of the residual error.

% W is FxP.  V is rxP, U is Fxr.  UV approximates W.  
% ENTERING_SHUM = 1
MAXSTEPS = num_it;
stable = 1;
Wdim = size(W);
F = Wdim(1);
P = Wdim(2);

if (sum(sum(INC)) < r*(F + P - r) + P)
% Return -2 is there's not enough information for this iterative method to work.
%   See Shum et al, for an explanation of this formula.
  error=-2;
  Wapp = [];
else
if INIT == 0
  V = rand(P,r);
else
  [Ui,Si,Vi] = svd(INIT);
  V = (Si(1:r,:)*Vi')';
end

Wm = W.*INC;

old_error = 999999999.0;
error = 77777777.0;
j = 1;

while (j <= MAXSTEPS) & (abs(error - old_error) > 1e-10)
  j = j+1;

  [U,s] = whole_invert(V, W, INC, r);
  stable = s & stable;
  [V,s] = whole_invert(U, W', INC', r);
  stable = s & stable;

  Wapp = U*V';
  old_error = error;
  error = sum(sum((Wm - Wapp.*INC).^2));

end

if j == MAXSTEPS+1
  [j, error, old_error];
  error = - error;
end

end  % if on conditioning of matrix.
