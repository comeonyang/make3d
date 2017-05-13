
function [At,c,K,yfeas] = createrandomdata(nF,nG,r,m)
% [At,c,K,yfeas] = createrandomdata(nF,nG,r,m);
%
% CREATERANDOMDATA randomly generates feasible rank constrained LMI problems 
% that can be used to test LMIRANK.
%
% Inputs:
%       nF  : the Fi matrices are size nF x nF        
%       nG  : the Gi matrices are size nG x nG
%       r   : rank bound
%       m   : number of varibales
%
% Outputs:
%       At      = -[vec(F1),...,vec(Fm); 
%                   vec(G1),...,vec(Gm)] 
%       c       = [vec(F0); 
%                  vec(G0)] 
%       K       : problem parameters (K.s and K.rank)
%       yfeas   : known feasible solution
%
% See also LMIRANK, LMIRANKTEST.

% Author Robert Orsi
% Feb 2005


%%%% Create y
y=randn(m,1); 
yfeas=y; % used for testing only

%%%% Create At
At=zeros(nF^2+nG^2,m);
c=zeros(nF^2+nG^2,1);
for i=1:m,
    Fi=randn(nF);
    Fi=sqrt(0.5)*(Fi+Fi');
    Gi=randn(nG);
    Gi=sqrt(0.5)*(Gi+Gi');
    At(:,i)=-[Fi(:); Gi(:)];
end

%% Create X=X'>=0
  X=randn(nF);
  X=X+X';
  [V,D]=eig(X);
  d=randn(nF,1);
  d=max(d,0);
  X=V*diag(d)*V';
  X=(X+X')/2; 
%% Create Y=Y'>=0, rank(Y)=r
  Y=randn(nG);
  Y=Y+Y';
  [V,D]=eig(Y);
  d=[rand(r,1);zeros(nG-r,1)];
  Y=V*diag(d)*V';
  Y=(Y+Y')/2;  
  
%%%% Create c  
c=At*y+[X(:);Y(:)];

%%%% Create K 
K.s=[nF nG];
K.rank=[nF r];
  
  
