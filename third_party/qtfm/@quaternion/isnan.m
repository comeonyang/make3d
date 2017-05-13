function tf = isnan(A)
% ISNAN  True for Not-a-Number.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(A)
    tf = isnan(x(A)) | isnan(y(A)) | isnan(z(A));
else
    tf = isnan(s(A)) | isnan(v(A));
end
