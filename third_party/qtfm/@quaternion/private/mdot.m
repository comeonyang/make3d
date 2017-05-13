function d = mdot(a, b)
% Matrix dot (scalar) product. Like dot, but matrix not elementwise.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))

if ~isa(a, 'quaternion') | ~isa(b, 'quaternion')
    error('mdot product is not defined for a quaternion and a non-quaternion.')
end

% This function is defined for full and pure quaternions, and combinations
% of full and pure, in which case we assume a zero scalar part for the pure
% argument.

if ispure(a) | ispure(b)
    % This covers the case where either or both are pure. We can ignore the
    % scalar part of the other, since it is implicitly multiplied by zero.
    
    d =               x(a) * x(b) + y(a) * y(b) + z(a) * z(b);
else
    d = s(a) * s(b) + x(a) * x(b) + y(a) * y(b) + z(a) * z(b);
end
