function b = repmat(a, m, n)
% REPMAT Replicate and tile an array.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 3, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 2
    if ispure(a)
        b = quaternion(repmat(x(a), m), ...
                       repmat(y(a), m), ...
                       repmat(z(a), m));
    else
        b = quaternion(repmat(s(a), m), ...
                       repmat(v(a), m));
    end
else
    if ispure(a)
        b = quaternion(repmat(x(a), m, n), ...
                       repmat(y(a), m, n), ...
                       repmat(z(a), m, n));
    else
        b = quaternion(repmat(s(a), m, n), ...
                       repmat(v(a), m, n));
    end
end
