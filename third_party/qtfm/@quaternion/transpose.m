function t = transpose(a)
% .'  Transpose.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout)) 

if ispure(a)
    t = quaternion(                 transpose(x(a)), transpose(y(a)), transpose(z(a)));
else
    t = quaternion(transpose(s(a)), transpose(x(a)), transpose(y(a)), transpose(z(a)));
end
