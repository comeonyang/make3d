function z = isocomplex(q)
% Construct a complex number from a quaternion, with the same modulus and
% argument as those of the quaternion.

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

% First extract the scalar part and the modulus of the vector part of the
% quaternion array q.

s =     scalar(q);
v = abs(vector(q));

% Either of both of s and v may be complex, in which case we cannot create
% an isomorphic complex number.

if isreal(s) & isreal(v)
  z = complex(s, v);
else
    error('Private function isocomplex called with complex quaternion argument')
end;

% Notes: the result will always have an argument in the range 0 to pi.
