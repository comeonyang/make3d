
function [At,c,K,yfeas,y,info] = lmiranktest
% [At,c,K,yfeas,y,info] = lmiranktest;
%
% LMIRANKTEST runs a test problem for LMIRANK. It uses randomly generated 
% feasible rank constrained LMI problems created by CREATERANDOMDATA.
%
% Problem parameters:
%       nLP     : no. of LP inequalities. This can be 0.
%       nF      : the Fi matrices are size nF x nF        
%       nG      : the Gi matrices are size nG x nG
%       r       : rank bound
%       m       : number of variables
%       maxiter : max. no. of iterations
%
% Outputs:  (LP ineq. data at the start of At and c is not shown)
%       At          = -[vec(F1),...,vec(Fm); 
%                       vec(G1),...,vec(Gm)] 
%       c           = [vec(F0); 
%                      vec(G0)] 
%       yfeas       : known feasible solution
%       y           : calculated solution
%       info.solved : 1 if a solution was found, 0 otherwise
%       info.cpusec : solution time
%       info.iters  : no. of iterations required to find a solution
%       info.gap    : constraint gap
%       info.rank   : ranks (with respect to tolerance pars.eps)
%
% See also LMIRANK, CREATERANDOMDATA.

% Author Robert Orsi
% Feb 2005


%%%% Random problem parameters
nLP=3;
nF=10;
nG=10;
r=5; 
m=20;  

%%%% Max. no. of iterations
maxiter=100;

%%%% Create random data
[At,c,K,yfeas] = createrandomdata(nF,nG,r,m);
%%%% Add LP ineq. constraints
At=[rand(nLP,m); At];
c=[(At(1:nLP,:)*yfeas+max(randn(nLP,1),0)); c];
K.l=nLP;

%%%% Call LMIRank
pars.maxiter=maxiter;
[y,info] = lmirank(At,c,K,pars);

