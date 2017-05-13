function d = tril(v, k)
% TRIL Extract lower triangular part.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout))

if nargin == 1
    if ispure(v)
        d = quaternion(            tril(x(v)), tril(y(v)), tril(z(v)));
    else
        d = quaternion(tril(s(v)), tril(x(v)), tril(y(v)), tril(z(v)));
    end
else
    if ispure(v)
        d = quaternion(               tril(x(v), k), tril(y(v), k), tril(z(v), k));
    else
        d = quaternion(tril(s(v), k), tril(x(v), k), tril(y(v), k), tril(z(v), k));
    end
end
