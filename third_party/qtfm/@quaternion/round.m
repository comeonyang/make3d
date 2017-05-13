function a = round(q)
% ROUND Round towards nearest integer.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout)) 

if ispure(q)
    a = quaternion(             round(x(q)), round(y(q)), round(z(q)));
else
    a = quaternion(round(s(q)), round(x(q)), round(y(q)), round(z(q)));
end
