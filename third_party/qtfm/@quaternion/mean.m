function m = mean(X, dim)
% MEAN   Average or mean value.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(X)
        m = quaternion(mean(x(X)), ...
                       mean(y(X)), ...
                       mean(z(X)));
    else
        m = quaternion(mean(s(X)), mean(v(X)));
    end
else
    if ispure(X)
        m = quaternion(mean(x(X), dim), ...
                       mean(y(X), dim), ...
                       mean(z(X), dim));
    else
        m = quaternion(mean(s(X), dim), mean(v(X), dim));
    end
end
