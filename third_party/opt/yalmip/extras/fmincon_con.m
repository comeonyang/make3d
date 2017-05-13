function [g,geq,dg,dgeq] = fmincon_con(x,params)

% Early bail for linear problems
if params.linearconstraints%isempty(params.interfacedata.evalMap) & isempty(params.interfacedata.Anonlinineq) & isempty(params.interfacedata.Anonlineq) & isequal( params.interfacedata.K.q,0) & isequal( params.interfacedata.K.s,0)
    g = [];
    geq = [];
    return
end

xevaled = zeros(1,length(params.interfacedata.c));
xevaled(params.linearindicies) = x;

% Experimental support for arbitrary functions
% nonlinear expressions inside sin() exp() etc
if ~isempty(params.interfacedata.evalMap)
    pp = params.monomtable(params.nonlinearindicies,:);
    xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^pp,2);

    for i = 1:length(params.interfacedata.evalMap)
        arguments = {params.interfacedata.evalMap{i}.fcn,xevaled(params.interfacedata.evalMap{i}.variableIndex)};
        arguments = {arguments{:},params.interfacedata.evalMap{i}.arg{2:end-1}};
        xevaled(params.interfacedata.evalVariables(i)) = feval(arguments{:});        
    end
end
pp = params.monomtable(params.nonlinearindicies,:);
xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^pp,2);

if params.nonlinearinequalities 
    g = params.interfacedata.Anonlinineq*xevaled(:)-params.interfacedata.bnonlinineq;
else
    g = [];
end

if params.nonlinearequalities
    geq = params.interfacedata.Anonlineq*xevaled(:)-params.interfacedata.bnonlineq;
else
    geq = [];
end

K = params.interfacedata.K;
top = 1;
if K.q(1) > 0
    for i = 1:length(K.q)
        Axcd = params.interfacedata.F_struc(top:top+K.q(i)-1,:)*[1;xevaled(:)];
        g = [g;-(Axcd(1)^2-norm(Axcd(2:end),2)^2)];
        top = top + K.q(i);
    end
end
if K.s(1) > 0
    for i = 1:length(K.s)
        CminusA = params.interfacedata.F_struc(top:top+K.s(i)^2-1,:)*[1;xevaled(:)];
        CminusA = reshape(CminusA,K.s(i),K.s(i));
        [R,p] = chol(CminusA);
        if p
            g = [g;-min(eig(CminusA))];
        else
            g = [g;-log(det(CminusA))];
        end
        top = top + K.s(i)^2;
    end
end

% 
% dg = [];
% dgeq = [];
% if params.SimpleNonlinearConstraints    
%     dg = [];    
%     allA = [params.interfacedata.Anonlineq;params.interfacedata.Anonlinineq];
%     for i = 1:length(params.linearindicies)
%         xevaled = zeros(1,length(params.interfacedata.c));
%         xevaled(params.linearindicies) = x;
%         mt = params.monomtable;
%         oldpower = mt(:,params.linearindicies(i));
%         mt(:,params.linearindicies(i)) = mt(:,params.linearindicies(i))-1;
%         xevaled = prod(repmat(xevaled,size(mt,1),1).^mt,2);        
%         xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
%         dg = [dg allA*xevaled'];
%     end     
%     dgeq = dg(1:size(params.interfacedata.Anonlineq,1),:)';
%     dg = dg(size(params.interfacedata.Anonlineq,1)+1:end,:)';
%     %full(dgeq')
% end

% 
% dg = [];
% dgeq = [];
% if params.SimpleNonlinearConstraints    
%     dg = [];    
%     allA = [params.interfacedata.Anonlineq;params.interfacedata.Anonlinineq];
%      mt = params.monomtable;
%     xe = zeros(1,length(params.interfacedata.c));
%     xe(params.linearindicies) = x;    
%     xx = xe;
%     xe = repmat(xe,size(mt,1),1).^mt;       
%    
%     for i = 1:length(params.linearindicies)
%        % xevaled = zeros(1,length(params.interfacedata.c));
%        % xevaled(params.linearindicies) = x;
% 
%         oldpower = mt(:,params.linearindicies(i));
%         newpower = oldpower-1;
% %        mt(:,params.linearindicies(i)) = mt(:,params.linearindicies(i))-1;
%         xevaled = xe;xevaled(params.linearindicies(i),:) = xx(params.linearindicies(i),:)
%         xevaledparams.linearindicies(i)) = xevaled(params.linearindicies(i),:).^newpower';
%         xevaled = prod(xevaled,2);        
%         xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
%         dg = [dg allA*xevaled'];
%     end     
%     dgeq = dg(1:size(params.interfacedata.Anonlineq,1),:)';
%     dg = dg(size(params.interfacedata.Anonlineq,1)+1:end,:)';
%     full(dgeq)
% end
% 



