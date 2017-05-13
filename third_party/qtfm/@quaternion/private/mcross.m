function c = mcross(a, b)
% Matrix cross (vector) product of two pure quaternions. Like cross but matrix
% not elementwise.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))

if ~isa(a, 'quaternion') | ~isa(b, 'quaternion')
    error('Matrix cross product is not defined for a quaternion and a non-quaternion.')
end

if ~ispure(a) | ~ispure(b)
    error('Mcross product is defined only for pure quaternions.')
end

c = quaternion(y(a) * z(b) - z(a) * y(b),   ...
               z(a) * x(b) - x(a) * z(b),   ...
               x(a) * y(b) - y(a) * x(b));
