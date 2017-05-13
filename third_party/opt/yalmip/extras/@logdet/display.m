function display(X)
%display           Overloaded

% Author Johan Löfberg 
% $Id: display.m,v 1.1 2004/06/17 08:40:08 johanl Exp $  


P = X.P;
classification = 'Logdet-term ';
[n,m] = size(P);
disp([classification num2str(n) 'x' num2str(m)]);
