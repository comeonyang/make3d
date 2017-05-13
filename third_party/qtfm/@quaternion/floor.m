function a = floor(q)
% FLOOR  Round towards minus infinity.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout)) 

if ispure(q)
    a = quaternion(             floor(x(q)), floor(y(q)), floor(z(q)));
else
    a = quaternion(floor(s(q)), floor(x(q)), floor(y(q)), floor(z(q)));
end
