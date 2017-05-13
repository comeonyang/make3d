function tf = isfinite(A)
% ISFINITE  True for finite elements.  
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(A)
    tf = isfinite(x(A)) & isfinite(y(A)) & isfinite(z(A));
else
    tf = isfinite(s(A)) | isfinite(v(A));
end
