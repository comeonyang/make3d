function issym=issymmetric(X)
%ISSYMMETRIC Check if variable is symmetric

% Author Johan Löfberg 
% $Id: issymmetric.m,v 1.6 2006/07/26 20:17:58 joloef Exp $  

n = X.dim(1);
m = X.dim(2);
if (n==m)
    % What are the linar indicies to the transposed matrices
    if isa(X.basis,'lazybasis')
        issym = 1;
    else
        ind = reshape(reshape(1:n^2,n,n)',n^2,1);
        issym = norm(X.basis-X.basis(ind,:),1)<1e-10;   
    end
else
    issym = 0;
end

