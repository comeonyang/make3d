function r = cat(dim, varargin)
% CAT Concatenate arrays.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(3, inf, nargin)), error(nargoutchk(0, 1, nargout))

a = varargin{1}; 
b = varargin{2};

if ispure(a) ~= ispure(b)
    error('Cannot concatenate mixture of pure and full quaternion arrays');
end

if ispure(a)
    r = quaternion(cat(dim, x(a), x(b)), ...
                   cat(dim, y(a), y(b)), ...
                   cat(dim, z(a), z(b)));
else
    r = quaternion(cat(dim, s(a), s(b)), ...
                   cat(dim, x(a), x(b)), ...
                   cat(dim, y(a), y(b)), ...
                   cat(dim, z(a), z(b)));
end

if nargin > 3
    r = cat(dim, r, varargin{3:end});    
end
