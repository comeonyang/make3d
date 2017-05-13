function varargout = recover(lmi_variables)
%RECOVER Create SDPVAR object using variable indicies

% Author Johan Löfberg
% $Id: recover.m,v 1.7 2005/02/28 15:49:44 johanl Exp $

if isempty(lmi_variables)
    varargout{1} = [];
else
    n = length(lmi_variables);
    i=1:n;
    if nargout <= 1
        varargout{1} = sdpvar(n,1,[],lmi_variables(:)',sparse(i,i+1,ones(n,1),n,n+1),0);
    else
        x = sdpvar(n,1,[],lmi_variables(:)',sparse(i,i+1,ones(n,1),n,n+1),0);
        for i = 1:length(lmi_variables)
            varargout{i} = x(i);
        end
    end
end