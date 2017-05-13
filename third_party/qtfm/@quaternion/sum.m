function t = sum(a, dim)
% SUM Sum of elements.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(a)
        t = quaternion(sum(x(a)), ...
                       sum(y(a)),....
                       sum(z(a)));
    else
        t = quaternion(sum(s(a)), sum(v(a)));
    end
else
    if dim == 'double' | dim == 'native'
        error('Parameters ''double'' or ''native'' are not implemented.');
    else
        if ispure(a)
            t = quaternion(sum(x(a), dim), ...
                           sum(y(a), dim), ...
                           sum(z(a), dim));
        else
            t = quaternion(sum(s(a), dim), sum(v(a), dim));
        end
    end
end
