function d = triu(v, k)
% TRIU Extract upper triangular part.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(v)
        d = quaternion(            triu(x(v)), triu(y(v)), triu(z(v)));
    else
        d = quaternion(triu(s(v)), triu(x(v)), triu(y(v)), triu(z(v)));
    end
else
    if ispure(v)
        d = quaternion(               triu(x(v), k), triu(y(v), k), triu(z(v), k));
    else
        d = quaternion(triu(s(v), k), triu(x(v), k), triu(y(v), k), triu(z(v), k));
    end
end
