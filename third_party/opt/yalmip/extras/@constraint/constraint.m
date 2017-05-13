function F = constraint(X,quantifier,Y)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: constraint.m,v 1.5 2005/06/02 13:40:00 joloef Exp $

superiorto('sdpvar');
superiorto('double');

if isa(X,'blkvar')
    X = sdpvar(X);
end
if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

% Try to evaluate
switch quantifier
case {'>','>='}
    Z = X - Y;
case {'<','<=','=='}
    Z = Y - X;
otherwise
    error('Quantifier not supported')
end

F.List={X,quantifier,Y};
F.Evaluated{1} = Z;

switch quantifier
case {'>','<'}
    F.strict(1) = 1;
case {'>=','<=','=='}
    F.strict(1) = 0;
otherwise
    error('Quantifier not supported')
end

F = class(F,'constraint');
	