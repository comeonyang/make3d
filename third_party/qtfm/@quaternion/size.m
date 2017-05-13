function [r, c] = size(q, dim)
% SIZE   Size of matrix.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 2, nargout))

xq = x(q); [xr, xc] = size(xq);
if ispure(q)
    if [xr, xc] ~= size(y(q)) | [xr, xc] ~= size(z(q))
        error('Sizes within object differ.')
    end
else
    if size(s(q)) ~= [xr, xc] | [xr, xc] ~= size(y(q)) | [xr, xc] ~= size(z(q))
        error('Sizes within object differ.')
    end
end

% In what follows, we use the size of the x component of the quaternion, not the size
% of the scalar part, since this could be empty. Otherwise the choice is arbitrary,
% since the code above checks that all components of q have the same size.

switch nargout
case 0
    switch nargin
    case 1
        size(xq)
    case 2
        size(xq, dim)
    end
case 1
    switch nargin
    case 1
        r = [xr, xc];
    case 2
        r = size(xq, dim);
    end
case 2
   switch nargin
   case 1
       [r, c] = size(xq);
   case 2
       [r, c] = size(xq, dim);
   end
otherwise
   error('Unhandled case.')
end



