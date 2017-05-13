function [F,properties,arguments] = mpower_internal(X,method,options,extstruct)
switch method
    case {'graph','milp'}
        F = set(extstruct.arg{1} == extstruct.var);
        arguments = extstruct.arg{1};
        properties = struct('convexity','none','monotonicity','none','definiteness','none');
    otherwise
        F = [];
        return
end
