function b = subsref(a, ss)
% SUBSREF Subscripted reference.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if length(ss) ~= 1
    error('Only one level of subscripting is currently supported for quaternions.');
    % See the notes below under structure indexing.
end

switch ss.type
case '()'
    if length(ss) ~= 1
        error('Multiple levels of subscripting are not supported for quaternions.')
    end
    
    % To implement indexing, we separate the quaternion into components, and then
    % compose the quaternion from the indexed components.
    
    if ispure(a)
        xa = x(a); ya = y(a); za = z(a);
        b = quaternion(xa(ss.subs{:}), ya(ss.subs{:}), za(ss.subs{:}));
    else
        sa = s(a); xa = x(a); ya = y(a); za = z(a);
        b = quaternion(sa(ss.subs{:}), xa(ss.subs{:}), ya(ss.subs{:}), za(ss.subs{:}));
    end
case '{}'
    error('Cell array indexing is not valid for quaternions.')
case '.'
    % Structure indexing.
    %
    % See some notes on this subject in the file subsasgn.m. Here, we would
    % have to support two levels of indexing, such as q.x(1,2) and q(1,2).x
    % which would give the same result. The Matlab help on subsref explains
    % how this would work. We could accept s, x, y, or z as field names and
    % return x(a) etc, but we would really need to pass a second level of
    % subscripting back to subsref recursively. The code below is an
    % interim step to support one level of structure indexing.
    %
    % 15-Sep-2005 Code contributed by T.A. Ell.
    %
    switch ss.subs
        case {'vector', 'v'}    
            b = vector(a); % maybe private function access would be better?
        case {'scalar', 's'}
            b = scalar(a);
        case {'x', 'I'}
            b = x(a);
        case {'y', 'J'}
            b = y(a);
        case {'z', 'K'}
            b = z(a);
        case 'imag'
            b = imag(a);
        case {'real'}
            b = real(a);
        otherwise
            error( ['Structure ''.', ss.subs, ''' is not a valid index']);
    end
otherwise
    error('subsref received an invalid subscripting type.')
end
