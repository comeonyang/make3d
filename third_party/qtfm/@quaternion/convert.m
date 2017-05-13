function p = convert(q, t)
% CONVERT Convert elements of a quaternion to another type.
%
% This function converts a quaternion array into a quaternion array with
% components of a different data type (e.g. int8, double, single).

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))

if ~ischar(t)
    error('Second parameter must be a string.')
end

f = str2func(t); % Construct a function handle from t, which must denote
                 % a function on the current Matlab path, so if it does
                 % not, an error will occur here.
if ispure(q)
    p = quaternion(         f(x(q)), f(y(q)), f(z(q)));
else
    p = quaternion(f(s(q)), f(x(q)), f(y(q)), f(z(q)));
end
