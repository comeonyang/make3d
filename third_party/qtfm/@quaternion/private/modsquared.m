function n = modsquared(q)
% Modulus squared of a quaternion.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

% Implementation note. We do not use q .* conj(q) because this would give
% a quaternion result with a zero vector part. It is better not to compute
% the vector part.

if ispure(q)
    n =            x(q).^2 + y(q).^2 + z(q).^2;
else
    n = s(q).^ 2 + x(q).^2 + y(q).^2 + z(q).^2;
end
