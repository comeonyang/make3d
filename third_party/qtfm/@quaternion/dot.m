function d = dot(a, b)
% DOT  Quaternion dot product.
% The dot (scalar) product of two quaternions is the sum of the products of
% the (s, x, y, z) components of the two quaternions.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))

if ~isa(a, 'quaternion') | ~isa(b, 'quaternion')
    error('Dot product is not defined for a quaternion and a non-quaternion.')
end

% This function is defined for full and pure quaternions, and combinations
% of full and pure, in which case we assume a zero scalar part for the pure
% argument.

if ispure(a) | ispure(b)
    % This covers the case where either or both are pure. We can ignore the
    % scalar part of the other, since it is implicitly multiplied by zero.
    
    d =                x(a) .* x(b) + y(a) .* y(b) + z(a) .* z(b);
else
    d = s(a) .* s(b) + x(a) .* x(b) + y(a) .* y(b) + z(a) .* z(b);
end
