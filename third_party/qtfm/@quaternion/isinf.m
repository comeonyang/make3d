function tf = isinf(A)
% ISINF  True for infinite elements.  
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(A)
    tf = isinf(x(A)) | isinf(y(A)) | isinf(z(A));
else
    tf = isinf(s(A)) | isinf(v(A));
end
