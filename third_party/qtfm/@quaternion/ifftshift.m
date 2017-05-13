function Y = ifftshift(X,dim)
% IFFTSHIFT Inverse FFT shift.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(X)
        Y = quaternion(                 ifftshift(x(X)), ifftshift(y(X)), ifftshift(z(X)));
    else
        Y = quaternion(ifftshift(s(X)), ifftshift(x(X)), ifftshift(y(X)), ifftshift(z(X)));
    end
else
    if ispure(X)
        Y = quaternion(                      ifftshift(x(X), dim), ifftshift(y(X), dim), ifftshift(z(X), dim));
    else
        Y = quaternion(ifftshift(s(X), dim), ifftshift(x(X), dim), ifftshift(y(X), dim), ifftshift(z(X), dim));
    end
end
