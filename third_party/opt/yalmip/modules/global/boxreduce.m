% *************************************************************************
% Bound strengthening
% *************************************************************************
function [p,feasible,vol_reduction] = boxreduce(p,upper,lower,lpsolver,options);

if options.bmibnb.lpreduce

    vol_start    = prod(p.ub(p.branch_variables)-p.lb(p.branch_variables));
    diag_before  =  sum(p.ub(p.branch_variables)-p.lb(p.branch_variables));

    [pcut,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver);
    diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
    iterations = 0;
    while (diag_after/(1e-18+diag_before) < 0.75    ) & feasible & iterations<4
        [pcut,feasible,lower] = lpbmitighten(pcut,lower,upper,lpsolver);
        diag_before = diag_after;
        diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
        iterations = iterations + 1;
    end

    % Clean up...
    for i = 1:length(pcut.lb)
        if (pcut.lb(i)>pcut.ub(i)) & (pcut.lb-pcut.ub < 1e-3)
            pcut.lb(i)=pcut.ub(i);
            pcut = updatenonlinearbounds(pcut,i);
        end
    end
    p.lb = pcut.lb;
    p.ub = pcut.ub;

    % Metric = (V0/V)^(1/n)
    vol_reduction = 1;%max(0,min(1,(prod(p.ub(p.branch_variables)-p.lb(p.branch_variables))/(1e-7+vol_start))^(1/length(p.branch_variables))));
    p.lb(p.lb<-1e12) = -inf;
    p.ub(p.ub>1e12) = inf;
else
    vol_reduction = 1;
    feasible = 1;
end
