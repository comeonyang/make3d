function [E, stable] = extend_matrix(M,INC,subM, rows, nonrows, r)
% rows indicates which rows of M and INC were used to find a solution.
%   subM is a fit to just these rows.  nonrows indicate rows of M
%   that still need to be fit.
E(rows,:) = subM;
stable = 1;

if ~isempty(nonrows)

  [u,s,v] = svd(subM);
  vp = v';
  basis = vp(1:r,:);
  for i = nonrows
     INCrow = INC(i,:);
     Mrow = M(i,:).*INCrow;
     INCmat = repeat(INCrow',r)';
     basisrows = basis.*INCmat;
     if rank(basisrows) ~= r 
       stable = 0

% putting zeros in for the elements that cannot be determined.
% to avoid dividing by a matrix that is not full rank.
       E(i,:) = Mrow;

%instead of zeros, choose a random value between -1 and 1 for each
%unfindable missing point. 
       missingpoints = find(INCrow == 0);
       E(i,missingpoints) = rand(1,length(missingpoints)).*2 - 1;
     else
       E(i,:) = (Mrow/basisrows)*basis;
     end
  end

end   
