function c = cross(a, b)
% CROSS (vector) product of two pure quaternions.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if ~isa(a, 'quaternion') | ~isa(b, 'quaternion')
    error('Cross product is not defined for a quaternion and a non-quaternion.')
end

if ~ispure(a) | ~ispure(b)
    error('Cross product is defined only for pure quaternions.')
end

c = quaternion(y(a) .* z(b) - z(a) .* y(b),   ...
               z(a) .* x(b) - x(a) .* z(b),   ...
               x(a) .* y(b) - y(a) .* x(b));
