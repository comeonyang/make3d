function [f,df] = fmincon_fun(x,params)

xevaled = zeros(1,length(params.interfacedata.c));
xevaled(params.linearindicies) = x;

% Experimental support for arbitrary functions
if ~isempty(params.interfacedata.evalMap)
    pp = params.monomtable(params.nonlinearindicies,:);
    xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^pp,2);
    for i = 1:length(params.interfacedata.evalMap)
        arguments = {params.interfacedata.evalMap{i}.fcn,xevaled(params.interfacedata.evalMap{i}.variableIndex)};
        arguments = {arguments{:},params.interfacedata.evalMap{i}.arg{2:end-1}};
        xevaled(params.interfacedata.evalVariables(i)) = feval(arguments{:});
    end
end

xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^params.monomtable(params.nonlinearindicies,:),2);

xevaled = xevaled(:);
f = params.interfacedata.c'*xevaled+xevaled'*params.interfacedata.Q*xevaled;

if params.SimpleLinearObjective
    df = params.interfacedata.c(params.linearindicies);
elseif params.SimpleQuadraticObjective
    df = params.interfacedata.c(params.linearindicies) + 2*params.interfacedata.Q(params.linearindicies,params.linearindicies)*x;
elseif params.SimpleNonlinearObjective
    df = [];
    for i = 1:length(params.linearindicies)
        xevaled = zeros(1,length(params.interfacedata.c));
        xevaled(params.linearindicies) = x;
        mt = params.monomtable;
        oldpower = mt(:,params.linearindicies(i));
        mt(:,params.linearindicies(i)) = mt(:,params.linearindicies(i))-1;
        xevaled = prod(repmat(xevaled,size(mt,1),1).^mt,2);
        xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
        df = [df;params.interfacedata.c'*xevaled'];
    end
    df = df + 2*params.interfacedata.Q(params.linearindicies,params.linearindicies)*x;
else
    df = [];
end
