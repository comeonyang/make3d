function  out = isequal(X,Y)
%ISEQUAL (overloaded)

% Author Johan Löfberg 
% $Id: isequal.m,v 1.2 2004/07/01 11:17:11 johanl Exp $   

if (isa(X,'sdpvar') & isa(Y,'sdpvar'))
    out = isequal(struct(X),struct(Y));
else
	out = 0;
end
	