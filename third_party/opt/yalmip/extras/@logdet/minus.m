function Z = minus(cx,P)
%display           Overloaded

% Author Johan Löfberg 
% $Id: minus.m,v 1.2 2004/09/19 21:33:49 johanl Exp $  

% Standard case c't-logdet(P)

if isa(P,'logdet') % sdpvr - logdet

    if prod(size(cx))>1
        error('Only scalar terms can be added to a logdet term');
    end
    
    if isa(cx,'logdet')
        error('Logdet objects can only be added');
    end
    
    Z = P;
    if isempty(P.cx)
        Z.cx = cx;
    else
        Z.cx = cx-P.cx;
    end
    Z.gain = -Z.gain;
else % logdet - cx
    temp = cx;
    cx = P;
    P = temp;
    
    Z = P;
    if isempty(P.cx)
        Z.cx = -cx;
    else
        Z.cx = P.cx-cx;
    end
end