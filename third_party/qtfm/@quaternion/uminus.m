function u = uminus(a)
% -  Unary minus.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout)) 

if ispure(a)
    u = quaternion(       -x(a), -y(a), -z(a));
else
    u = quaternion(-s(a), -x(a), -y(a), -z(a));
end
