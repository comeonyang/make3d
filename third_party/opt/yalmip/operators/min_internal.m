function [F,properties,arguments] = min_internal(X,method,options,extstruct)
switch method
    case 'graph'
        arguments=[];
        F = set([]);
        for j = 1:length(extstruct.arg)
            F = F + set(extstruct.arg{j} - extstruct.var);
            arguments= [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none');
    case 'milp'
        arguments = [];
        F = set([]);
        t = extstruct.var;
        for j = 1:length(extstruct.arg) % MAX(x,y)
            X = extstruct.arg{j};
            X = reshape(X,length(X),1);
            [M,m] = derivebounds(X);
            n = length(X);
            d = binvar(n,1);
            F = F + set(sum(d)==1);
            F = F + set(-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));

            kk = [];
            ii = [];
            for i = 1:n
                k = [1:1:i-1 i+1:1:n]';
                ii = [ii;repmat(i,n-1,1)];
                kk = [kk;k];
                Mm = M(k)-m(i);
            end
            xii = extsubsref(X,ii);
            dii = extsubsref(d,ii);
            xkk = extsubsref(X,kk);
            F = F + set(xii <= xkk+(M(ii)-m(kk)).*(1-dii));
            arguments = [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','exact','monotonicity','exact','definiteness','none');

    otherwise
        F = [];
        return
end