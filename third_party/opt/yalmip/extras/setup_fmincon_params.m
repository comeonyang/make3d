function constant_data = setup_fmincon_params(interfacedata)

monomtable = interfacedata.monomtable;
nonlinearindicies = interfacedata.nonlinearindicies;
linearindicies = interfacedata.linearindicies;

constant_data.interfacedata = interfacedata;
constant_data.monomtable = monomtable;
constant_data.nonlinearindicies = nonlinearindicies;
constant_data.linearindicies = linearindicies;

% Figure out if YALMIP easily can compute the gradient of the objective
% This will done completely general later
constant_data.SimpleLinearObjective = 0;
constant_data.SimpleQuadraticObjective = 0;
constant_data.SimpleNonlinearObjective = 1;
constant_data.SimpleNonlinearConstraints = 0;
if isempty(interfacedata.evalMap)
    if nnz(interfacedata.c(nonlinearindicies)) == 0
        if (nnz(interfacedata.Q)==0)
            constant_data.SimpleLinearObjective = 1;
        else
            if nnz(interfacedata.Q(nonlinearindicies,nonlinearindicies))==0
                constant_data.SimpleQuadraticObjective = 1;
            end
        end
    end
    if isequal(interfacedata.K.s,0) & isequal(interfacedata.K.q,0) & isequal(interfacedata.K.r,0)
        constant_data.SimpleNonlinearConstraints = 1;
    end
else
    constant_data.SimpleNonlinearObjective = 0;
end
constant_data.nonlinearineval = 0;
% Check if there are any nonlinear expression in evaluation based operators
if ~isempty(interfacedata.evalMap)
    temp = [];
    for i = 1:length(interfacedata.evalMap)
        temp = [temp interfacedata.evalMap{i}.variableIndex];
    end
    if any(interfacedata.variabletype(temp))
        constant_data.nonlinearineval = 1;
    end
end


