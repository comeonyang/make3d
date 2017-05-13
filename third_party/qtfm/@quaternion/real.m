function p = real(q)
% REAL   Real part of a quaternion.
% (Quaternion overloading of standard Matlab function.)
%
% This function returns the quaternion that is the real
% part of q.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(q)
    p = quaternion(real(x(q)), real(y(q)), real(z(q)));
else
    p = quaternion(real(s(q)), real(v(q)));
end

