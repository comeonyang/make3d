function a = ceil(q)
% CEIL   Round towards plus infinity.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout)) 

if ispure(q)
    a = quaternion(            ceil(x(q)), ceil(y(q)), ceil(z(q)));
else
    a = quaternion(ceil(s(q)), ceil(x(q)), ceil(y(q)), ceil(z(q)));
end
