function q = mtimes(l, r)
% *   Matrix multiply.
% (Quaternion overloading of standard Matlab function.)
 
% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if isa(l, 'quaternion') & isa(r, 'quaternion')
    
    % Both arguments are quaternions. There are four cases to handle,
    % dependent on whether the arguments are pure or full.
    
    if ispure(l) & ispure(r)
        q = -mdot(l, r) + mcross(l, r);
    elseif ispure(l)
        q =                 l  * s(r) +             +   l  * v(r);
    elseif ispure(r)
        q =                             s(l) *   r  + v(l) *   r;
    else
        q = s(l) * s(r) + v(l) * s(r) + s(l) * v(r) + v(l) * v(r);
    end
   
else
    % One of the arguments is not a quaternion. If it is numeric, we can handle it,
    % and we must because the code above requires us to multiply scalar parts by
    % vector parts using a recursive call to this function.
    
    if isa(l, 'quaternion') & isa(r, 'numeric')
        if ispure(l)
            q = quaternion(           x(l) * r, y(l) * r, z(l) * r);
        else
            q = quaternion(s(l) * r, x(l) * r, y(l) * r, z(l) * r);
        end
    elseif isa(l, 'numeric') & isa(r, 'quaternion')
        if ispure(r)
            q = quaternion(           l * x(r), l * y(r), l * z(r));
        else
            q = quaternion(l * s(r), l * x(r), l * y(r), l * z(r));
        end
    else
        error('Matrix multiplication of a quaternion by a non-numeric is not implemented.')
    end
end



