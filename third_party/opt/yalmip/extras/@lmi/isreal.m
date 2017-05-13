function YESNO = isreal(F)
%ISREAL (overloaded)

% Author Johan Löfberg 
% $Id: isreal.m,v 1.2 2005/02/04 10:10:27 johanl Exp $   

if isempty(F.clauses)
    YESNO = 1;
else
    YESNO = 1;
    i = 1;
    while (i<=length(F.clauses)) & YESNO
        Fi = F.clauses{i};
        YESNO =  YESNO & is(Fi.data,'real');
        i  =i+1;
    end   
end

