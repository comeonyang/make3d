function Z = plus(X,Y)
%display           Overloaded

% Author Johan Löfberg 
% $Id: plus.m,v 1.2 2004/09/19 21:33:49 johanl Exp $  

% LOGDET + SDPVAR
if isa(Y,'logdet')
    Z = X;
    X = Y;
    Y = Z;
end

if prod(size(Y))>1
    error('Only scalar terms can be added to a logdet term');
end

if isa(Y,'logdet')
    Z = X;
    Z.P = blkdiag(Z.P,Y.P);
    return
end


Z = X;
if isempty(Z.cx)
    Z.cx = Y;
else
    Z.cx = plus(Z.cx,Y);
end