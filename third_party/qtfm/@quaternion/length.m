function n = length(q)
% LENGTH   Length of vector.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This was written because length appears not to work correctly for arrays
% of structs, as can be demonstrated by:
%
% a.x = randn(5,1), length(a)
%
% which gives 1, and
%
% a = randn(5,1), length(a)
%
% which gives 5, as expected.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(q)
    if length(x(q)) ~= length(y(q)) | length(x(q)) ~= length(z(q))
        error('Sizes within object differ.')
    end
else
    if length(s(q)) ~= length(x(q)) | length(s(q)) ~= length(y(q)) | length(s(q)) ~= length(z(q))
        error('Sizes within object differ.')
    end
end

% In what follows, we use the length of the x component, not the length of the scalar part,
% since this could be empty. Otherwise the choice is arbitrary, since the code above checks
% that all components of q have the same length.

n = length(x(q));
