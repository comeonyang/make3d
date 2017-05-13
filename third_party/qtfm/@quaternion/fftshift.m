function Y = fftshift(X,dim)
% FFTSHIFT Shift zero-frequency component to center of spectrum.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if nargin == 1
    if ispure(X)
        Y = quaternion(                fftshift(x(X)), fftshift(y(X)), fftshift(z(X)));
    else
        Y = quaternion(fftshift(s(X)), fftshift(x(X)), fftshift(y(X)), fftshift(z(X)));
    end
else
    if ispure(X)
        Y = quaternion(                     fftshift(x(X), dim), fftshift(y(X), dim), fftshift(z(X), dim));
    else
        Y = quaternion(fftshift(s(X), dim), fftshift(x(X), dim), fftshift(y(X), dim), fftshift(z(X), dim));
    end
end
