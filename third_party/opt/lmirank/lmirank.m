
function [y,info] = lmirank(At,c,K,pars,y0)
% [y,info] = lmirank(At,c,K,pars,y0);
%
% LMIRANK can be used to try to solve rank constrained LMI problems such as
%
%       F(y)>=0,                                (1)
%       G(y)>=0,  rank G(y)<=r.                 (2)
%
% More precisely it can be used to try to solve feasibility problems
% involving any number of LMI constraints and one or more rank constraints.
%
% LMI data is entered in standard SeDuMi format. Rank constraints are entered
% using K.rank. For example, 
%
%   K.s=[4 7 6]; K.rank=[4 4 6];
%
% specifies 3 LMI constraints with the 2nd LMI constrained to have rank <= 4.
%
% LP inequality constraints can also be included.
%
% LMIRANK can be called using any of following
%    [y,info] = lmirank(At,c,K,pars,y0);
%    [y,info] = lmirank(At,c,K);
%    [y,info] = lmirank(At,c,K,pars);
%    [y,info] = lmirank(At,c,K,y0);
%
% Inputs:
%       At,c,K.l,K.s    : data in SeDuMi format 
%       K.rank          : rank constraints 
% Optional Inputs:
%       pars.maxiter    : max. no. of iterations, default is 1000
%       pars.eps        : constraint tolerance, default is 1e-9
%       pars.fid        : set to 0 to suppress on-screen output. The default 
%                         is 1 which displays on-screen output.
%       pars.itermod    : output results to screen every pars.itermod 
%                         iterations, default is 1    
%       y0              : initial condition 
%                         (If y0 is not given, the trace heuristic is
%                          used to initialize the algorithm. SeDuMi is
%                          used to do this calculation and hence must
%                          be installed if y0 is not given.) 
% Outputs:
%       y           : solution 
%       info.solved : 1 if a solution was found, 0 otherwise
%       info.cpusec : solution time
%       info.iters  : no. of iterations required to find a solution
%       info.gap    : constraint gap
%       info.rank   : ranks (with respect to tolerance pars.eps)
%
% The algorithm is described in 
%       R. Orsi, U. Helmke, and J. B. Moore. A Newton-like method for solving 
%       rank constrained linear matrix inequalities. In Proceedings of the 
%       43rd IEEE Conference on Decision and Control, pages 3138-3144, 
%       Paradise Island, Bahamas, 2004.  
%
% Algorithm notes 
%   1. Projection onto the set B (see Section V.C of the paper) is now done 
%      in a much simpler manner. (It still results in a linearly constrained 
%      least squares problem.)
%   2. The above paper describes only the basic case given by (1) and (2) 
%      above, i.e., it does not deal with multiple non-rank constrained LMIs,
%      multiple rank constrained LMIs, or LP inequality constraints.
%
%Feedback should be sent to robert.orsi@anu.edu.au

% Author Robert Orsi
% Feb 2005


t=cputime;

%%%% If no LP ineq. constraints are present, set K.l=0
if ~isfield(K,'l')
    K.l=0;
end

%%%% Set unspecified pars, calculate an initial condition if none given
if nargin==3
    pars.fid=1;
    y = trheuristic(At,c,K);
end
if nargin==4
    if isstruct(pars)
        y = trheuristic(At,c,K);   
    else
        y=pars;
        pars.fid=1;
    end
end
if nargin==5
    y=y0;
end
if ~isfield(pars,'maxiter') pars.maxiter=1000; end
if ~isfield(pars,'eps') pars.eps=1e-9; end
if ~isfield(pars,'fid') pars.fid=1; end
if ~isfield(pars,'itermod') pars.itermod=1; end 

%%%% Initialize X1
X1=c-At*y;
         
%%%% Calculate m
m=size(At,2);

%%%% Print output header
fprintf(pars.fid,'\nLMIRank by Robert Orsi, 2005.\n');
fprintf(pars.fid,'maxiter = %5d  |   rank bounds\n',pars.maxiter);
fprintf(pars.fid,'eps = %0.2e  | ',pars.eps);
fprintf(pars.fid,'%5d',K.rank);
fprintf(pars.fid,'\n iter :     gap      ranks\n');

%%%% MAIN LOOP
for iters=1:pars.maxiter,
    
    %% Calculate DX1        : eigs of X1 
    %%
    %%           X2LP       : proj. of LP ineq. component of X1
    %%           X2LPindex  : index of zeroes of X2LP 
    %%           DX2        : eigs of X2. 
    %%                      : Sorted, non-neg., rank constained version of DX1   
    %%           rankX2     : ranks of X2
    %%
    %%           Attrans    : At after transformation
    %%           ctrans     : c after transformation
    %%           Attrans2   : 
    %%           ctrans2    :
    DX1=zeros(sum(K.s),1);
    %
    X2LP=max(X1(1:K.l),0);
    X2LPindex=find(X2LP==0);
    DX2=zeros(sum(K.s),1);
    rankX2=zeros(length(K.s),1); 
    %
    Attrans=zeros(size(At));
    Attrans(1:K.l,:)=At(1:K.l,:);
    ctrans=zeros(size(c));
    ctrans(1:K.l)=c(1:K.l);
    Attrans2=zeros(size(At)); % incorrect size, truncated later
    Attrans2(1:length(X2LPindex),:)=Attrans(X2LPindex,:);
    ctrans2=zeros(size(c)); % incorrect size, truncated later
    ctrans2(1:length(X2LPindex))=ctrans(X2LPindex); 
    %
    index=0;
    index2=0;    
    index3=length(X2LPindex);
    for j=1:length(K.s),
        [V,D]=eig(reshape(X1(K.l+index2+1:K.l+index2+K.s(j)^2),K.s(j),K.s(j)));
        DX1(index+1:index+K.s(j))=diag(D);
        [Dsort,I]=sort(-diag(D));
        Dsort=-Dsort;
        V=V(:,I);
        for i=1:K.s(j),    
            if (rankX2(j)<K.rank(j)) & (Dsort(i)>0)
                rankX2(j)=rankX2(j)+1;
                DX2(index+i)=Dsort(i);
            end    
        end
        for k=1:m,
            Q=V'*reshape(At(K.l+index2+1:K.l+index2+K.s(j)^2,k),K.s(j),K.s(j))*V;
            Q=(Q+Q')/2;
            Attrans(K.l+index2+1:K.l+index2+K.s(j)^2,k)=Q(:); 
            Attrans2(index3+1:index3+(K.s(j)-rankX2(j))^2,k)=...
                reshape(Q(rankX2(j)+1:end,rankX2(j)+1:end),(K.s(j)-rankX2(j))^2,1);
        end
        Q=V'*reshape(c(K.l+index2+1:K.l+index2+K.s(j)^2),K.s(j),K.s(j))*V;
        Q=(Q+Q')/2;
        ctrans(K.l+index2+1:K.l+index2+K.s(j)^2)=Q(:); 
        ctrans2(index3+1:index3+(K.s(j)-rankX2(j))^2)=...
            reshape(Q(rankX2(j)+1:end,rankX2(j)+1:end),(K.s(j)-rankX2(j))^2,1);   
        %
        index=index+K.s(j);
        index2=index2+K.s(j)^2;
        index3=index3+(K.s(j)-rankX2(j))^2;
    end 
    %% Truncate Attrans2 and ctrans2 to correct sizes
    q=length(X2LPindex);
    for j=1:length(K.s),
        q=q+(K.s(j)-rankX2(j))^2;    
    end
    Attrans2=Attrans2(1:q,1:m);
    ctrans2=ctrans2(1:q,1);
    
    %% BREAK if a solution has been found
    %% Decision is based on X1, not X2
    bound= (K.l+length(DX1)) * max([abs(X1(1:K.l)); abs(DX1)]) * eps * 10;
    bound= max(bound,pars.eps);
    gap=min(DX1);
    if K.l~=0
        gap=min([gap; X1(1:K.l)]);
    end
    breakflag=(gap >= -bound);
    rankX1=zeros(1,length(K.s));
    index=0;
    for j=1:length(K.s),
        rankX1(j)=sum(abs(DX1(index+1:index+K.s(j)))>bound);
        breakflag=breakflag & (rankX1(j) <= K.rank(j));
        index=index+K.s(j); 
    end  
    if (pars.fid)&(mod(iters,pars.itermod)==0) %% Output iters, gap & ranks 
        fprintf(pars.fid,'%5d : %0.2e ',iters,gap);
        fprintf(pars.fid,'%5d',rankX1);
        fprintf(pars.fid,'\n');
    end
    if breakflag
        break
    end    
    
    %% Calculate SVD and rank of Attrans2
    [U,S,V]=svd(Attrans2);
    s=diag(S);
    tol=max(size(Attrans2))*s(1)*eps;
    rankAttrans2=sum(s>tol);
    
    %% Calculate y
    if rankAttrans2 == size(Attrans2,2),
        y=Attrans2\ctrans2; 
    else 
        if size(Attrans2,1)~=size(Attrans2,2)
            warning off;
            y0=Attrans2\ctrans2;
            warning on;
        else
            y0=pinv(Attrans2)*ctrans2;    
        end
        W=V(:,rankAttrans2+1:end);
        e=zeros(size(c));
        e(1:K.l)=X2LP;
        index=0;
        index2=K.l;
        for j=1:length(K.s),
            P=diag(DX2(index+1:index+K.s(j)));
            e(index2+1:index2+K.s(j)^2)=P(:);
            index=index+K.s(j);
            index2=index2+K.s(j)^2;
        end
        q=(Attrans*W)\( ctrans - Attrans*y0 - e );
        y=y0+W*q;
    end
      
    %% Calculate new X1
    X1=c-At*y;
     
end

info.solved=breakflag;
info.cpusec=cputime-t;
info.iters=iters;
info.gap=gap;           
info.rank=rankX1;

fprintf(pars.fid,'iters solved   seconds\n');
fprintf(pars.fid,'%5d %4d %11.1e\n',info.iters,info.solved,info.cpusec);
