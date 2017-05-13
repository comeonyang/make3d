function p = preprocess_eval_bounds(p);
if ~isempty(p.evalVariables)
    for i = 1:length(p.evalMap)
        arg = p.evalMap{i}.variableIndex;
        xL = p.lb(arg);
        xU = p.ub(arg);
        switch p.evalMap{i}.fcn
            case 'exp'
                p.lb(p.evalVariables(i)) = max([0  p.lb(p.evalVariables(i)) exp(xL)]);
                p.ub(p.evalVariables(i)) = min([p.ub(p.evalVariables(i)) exp(xU)]);
            case {'cos','sin'}
                if isequal(p.evalMap{i}.fcn,'cos')
                    xL = xL + pi/2;
                    xU = xU + pi/2;
                end

                neg = 0;
                if xL > 0
                    n = floor((xL/(2*pi)));
                    xL = xL - n*2*pi;
                    xU = xU - n*2*pi;
                else
                    n = floor((-xL/(2*pi)));
                    xL = xL + n*2*pi;
                    xU = xU + n*2*pi;
                    neg = 1;
                end
                yL = sin(xL);
                yU = sin(xU);
                L = min([yL yU]);
                U = max([yL yU]);
                if xL<pi/2 & xU>pi/2
                    U = 1;
                end
                if xL < 3*pi/2 & xU > 3*pi/2
                    L = -1;
                end
                if neg
                    t = L;
                    L = -U;
                    U = -t;
                end
                p.lb(p.evalVariables(i)) = max([p.lb(p.evalVariables(i)) L]);
                p.ub(p.evalVariables(i)) = min([p.ub(p.evalVariables(i)) U]);
                
            case 'tan'                
                n1 = fix((xL+pi/2)/(pi));
                n2 = fix((xU+pi/2)/(pi));
                if n1==n2
                    L = tan(xL);
                    U = tan(xU);
                else
                    L = -inf;
                    U = inf;                    
                end
                p.lb(p.evalVariables(i)) = max([p.lb(p.evalVariables(i)) L]);
                p.ub(p.evalVariables(i)) = min([p.ub(p.evalVariables(i)) U]);
                    

            case 'sdpfun'
                % We know nothing...
                % To get some kind of bounds, we jsut sample the function
                % and pick the min and max from there
                % this only works for simple functions
                if length(xL)>1
                    disp([p.evalMap{i}.fcn ' is not supported in the global solver (only scalar functions support)'])
                    error([p.evalMap{i}.fcn ' is not supported in the global solver'])
                end
                xtest = linspace(xL,xU,100);
                values = feval(p.evalMap{i}.arg{2},xtest);
                [minval,minpos] = min(values);
                [maxval,maxpos] = min(values);
                xtestmin = linspace(xtest(max([1 minpos-5])),xtest(min([100 minpos+5])),100);
                xtestmax = linspace(xtest(max([1 maxpos-5])),xtest(min([100 maxpos+5])),100);
                values1 = feval(p.evalMap{i}.arg{2},xtestmin);
                values2 = feval(p.evalMap{i}.arg{2},linspace(xL,xU,100));
                p.lb(p.evalVariables(i)) = max([p.lb(p.evalVariables(i)) min([values1 values2])]);
                p.ub(p.evalVariables(i)) = min([p.ub(p.evalVariables(i)) max([values1 values2])]);
            otherwise
                disp([p.evalMap{i}.fcn  ' is not supported in the global solver'])
                error([p.evalMap{i}.fcn ' is not supported in the global solver'])
        end
    end
end