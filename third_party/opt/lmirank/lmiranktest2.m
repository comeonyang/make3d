
function [A,Bperp,Ctperp,y,info,X,Y] = lmiranktest2
% [A,Bperp,Ctperp,y,info,X,Y] = lmiranktest2;
%
% LMIRANKTEST2 runs a test problem for LMIRANK. 
% The problem is a two-mass-spring, reduced order output feedback problem 
% considered in
%
%       R. Orsi, U. Helmke, and J. B. Moore. A Newton-like method for solving 
%       rank constrained linear matrix inequalities. In Proceedings of the 
%       43rd IEEE Conference on Decision and Control, pages 3138-3144, 
%       Paradise Island, Bahamas, 2004.  
%
% It is of the form:
% 
%       Find symmetric matrices X and Y such that
%
%        -Bperp*(A*X+X*A')*Bperp' > 0,          (1)
%      -Ctperp*(Y*A+A'*Y)*Ctperp' > 0,          (2)
%                      [X I; I Y] >= 0,         (3)
%                  rank[X I; I Y] <= n + k.     (4)
%% The code first uses YALMIP to convert constraints (1), (2) and (3) into 
% the SeDuMi data format. (This is the format used by LMIRANK.) After 
% specifying the rank constraint, LMIRANK is called. Using the solution of 
% LMIRANK, X and Y are reconstructed using YALMIP. Finally, the minimum 
% eigenvalues of the left hand sides of (1), (2) and (3) are displayed, 
% along with the rank in (4).
%
% See also LMIRANK.

% Author Robert Orsi
% Feb 2005


%%%% Problem data
alpha=0.2; 
A=[ 0  0 1 0; 0  0 0 1; -1  1 0 0; 1 -1 0 0];
A=A+alpha*eye(size(A));
Bperp=[1 0 0 0; 0 1 0 0; 0 0 0 1 ];
Ctperp=[1 0 0 0; 0 0 1 0; 0 0 0 1];
epsilon = 1e-6; 
k=2; % controller order

%%%% Determine various sizes
n=size(A,2);
m=size(Bperp,1);
p=size(Ctperp,1);

%%%% Create A,c,K using YALMIP
yalmip('clear');
X=sdpvar(n,n);
Y=sdpvar(n,n);
M1=-Bperp*(A*X+X*A')*Bperp';
M2=-Ctperp*(Y*A+A'*Y)*Ctperp';
M3=[X eye(n); eye(n) Y];
F=set( M1 > epsilon*eye(m) ) + set( M2 > epsilon*eye(p) ) + set( M3 > 0 );
[model,recoverymodel] = export(F,[],sdpsettings('solver','sedumi'));
At=model.A; 
c=model.C; 
K.s=model.K.s;

%%%% Specify rank constraint
K.rank=[m p (n+k)];

%%%% Call LMIRank
pars.itermod=10;   
[y,info] = lmirank(At,c,K,pars);

%%%% Reconstruct X and Y from y using YALMIP
setsdpvar(recover(recoverymodel.used_variables),y);
X=double(X);
Y=double(Y);

%%%% Testing 
M1=double(M1);
M2=double(M2);
M3=double(M3);
disp(' ');
disp(['min eig M1 =  ',num2str(min(eig(M1)),2)]);
disp(['min eig M2 =  ',num2str(min(eig(M2)),2)]);
disp(['min eig M3 =  ',num2str(min(eig(M3)),2)]);
disp(['rank M3 =  ',num2str(rank(M3),2)]);
disp(' ');

