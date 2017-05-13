function [p,feasible,vol_reduction] = domain_reduction(p,upper,lower,lpsolver);
% This is just too expensive
t1 = p.binary_variables;
t2 = p.integer_variables;
p.binary_variables = [];
p.integer_variables = [];
if ~p.options.bmibnb.lpreduce | (size(p.lpcuts,1)==0)
    vol_reduction = 1;
    p.feasible = 1;

    p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables));
    p.ub(p.integer_variables) = floor(p.ub(p.integer_variables));
    p.lb(p.binary_variables) = ceil(p.lb(p.binary_variables));
    p.ub(p.binary_variables) = floor(p.ub(p.binary_variables));

else
    [p,p.feasible,vol_reduction] =  boxreduce(p,upper,lower,lpsolver,p.options);
end
p.binary_variables  = t1;
p.integer_variables = t2;