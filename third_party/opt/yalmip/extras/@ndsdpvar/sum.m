function X = sum(varargin)
% SUM (overloaded)

% Author Johan Löfberg
% $Id: sum.m,v 1.4 2006/07/25 19:02:39 joloef Exp $

Y = varargin{1};
X = Y;
X.basis = [];
for i = 1:size(Y.basis,2)
    base = reshape(full(Y.basis(:,i)),X.dim);
    base = sum(base,varargin{2:end});
    X.basis = [X.basis sparse(base(:))];
end
X.dim = size(base);
X.conicinfo = [0 0];