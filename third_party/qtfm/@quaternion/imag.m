function p = imag(q)
% IMAG   Imaginary part of a quaternion.
% (Quaternion overloading of standard Matlab function.)
%
% This function returns the quaternion that is the imaginary
% part of q. If q is a real quaternion, it returns zero.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(q)
    p = quaternion(imag(x(q)), imag(y(q)), imag(z(q)));
else
    p = quaternion(imag(s(q)), imag(v(q)));
end

