function d = diag(v, k)
% DIAG Diagonal matrices and diagonals of a matrix.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(v)
        d = quaternion(            diag(x(v)), diag(y(v)), diag(z(v)));
    else
        d = quaternion(diag(s(v)), diag(x(v)), diag(y(v)), diag(z(v)));
    end
else
    if ispure(v)
        d = quaternion(               diag(x(v), k), diag(y(v), k), diag(z(v), k));
    else
        d = quaternion(diag(s(v), k), diag(x(v), k), diag(y(v), k), diag(z(v), k));
    end
end
