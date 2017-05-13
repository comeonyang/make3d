function [F,strict] = getlist(X)
% Internal class for constraint lists

% Author Johan Löfberg
% $Id: getlist.m,v 1.2 2005/05/01 11:52:01 joloef Exp $

superiorto('sdpvar');
superiorto('double');

F = X.Evaluated;
strict = X.strict;

% FIX : treat equalities better
if isequal(X.List{2},'==')
    F{1}=sethackflag(F{1},3);
end
