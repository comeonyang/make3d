function q = plus(l, r)
% +   Plus.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

% Three cases have to be handled:
%
% l is a quaternion, r is not,
% r is a quaternion, l is not,
% l and r are both quaternions.
%
% The complication here is that we can't call ispure() or s() etc
% for a parameter that is not a quaternion, so we have to handle
% these cases differently.

if isa(l, 'quaternion') & isa(r, 'quaternion')

  % The scalar part could be empty, and we have to handle this
  % specially, because [] + x is [] and not x, for some reason.

  if ispure(l) & ispure(r)
    q = quaternion(             x(l) + x(r), y(l) + y(r), z(l) + z(r));
  elseif ispure(l)
    q = quaternion(       s(r), x(l) + x(r), y(l) + y(r), z(l) + z(r));
  elseif ispure(r)
    q = quaternion(s(l),        x(l) + x(r), y(l) + y(r), z(l) + z(r));
  else
    q = quaternion(s(l) + s(r), x(l) + x(r), y(l) + y(r), z(l) + z(r));
  end

elseif isa(l, 'quaternion') & isa(r, 'numeric')

  if ispure(l)
    q = quaternion(       r, x(l), y(l), z(l));
  else
    q = quaternion(s(l) + r, x(l), y(l), z(l));
  end

elseif isa(l, 'numeric') & isa(r, 'quaternion')

  if ispure(r)
    q = quaternion(l       , x(r), y(r), z(r));
  else
    q = quaternion(l + s(r), x(r), y(r), z(r));
  end

else
  error('Unhandled parameter types in function +/plus')
end


