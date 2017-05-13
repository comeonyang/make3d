
function y = trheuristic(At,c,K)
% y = trheuristic(At,c,K);
%
% TRHEURISTIC uses SeDuMi to solve
%
%   min     tr(G(y))
%   s.t.    F(y) := F0 + y1*F1 + ... + ym*Fm >= 0,      
%           G(y) := G0 + y1*G1 + ... + ym*Gm >= 0.
%
% This is the trace heuristic for trying to solve rank constrained LMIs.
%
% TRHEURISTIC is used by LMIRANK for initialization.
%
% Inputs:
%       At      = -[vec(F1),...,vec(Fm); 
%                   vec(G1),...,vec(Gm)] 
%       c       = [vec(F0); 
%                  vec(G0)] 
%       K       : problem parameters (K.s and K.rank)
% 
% Outputs:
%       y       : solution   
%
% NOTE: It is possible to have multiple 'F' LMIs and multiple 'G' LMIs.
%       For simplicity, the description above assumes there is only 
%       one 'F' LMI and one 'G' LMI. If there are multiple 'G' LMIs, 
%       i.e., if K.s(j)~=K.rank(j) for more than one j, then only the 
%       last 'G' LMI is treated as such; the others are treated as 
%       'F' LMIs. LP inequality constraints can also be included.
%
% See also LMIRANK.

% Author Robert Orsi
% Feb 2005


%%%% SeDuMi:   
%   min     -b'*y
%   s.t.    c-At*y \in K 

%%%% Create b
m=size(At,2);
b=zeros(m,1);
if isfield(K,'l')
    index=K.l;
else    
    index=0;
end
for j=1:length(K.s),
    if K.s(j)~=K.rank(j)
        for i=1:m,
            P=reshape(At(index+1:index+K.s(j)^2,i),K.s(j),K.s(j));
            b(i)=trace(P);
        end
    end
    index=index+K.s(j)^2;
end

%%%% Set parameters
pars.fid=0; % Suppress SeDuMi on-screen output.
pars.eps=0; % was 1e-14.  
            % pars.eps=0 means SeDuMi runs as long as it can make progress.

%%%% Call SeDuMi          
[x,y,info] = sedumi(At',b,c,K,pars);




