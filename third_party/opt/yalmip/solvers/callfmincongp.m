function output = callfmincongp(interfacedata)

% Author Johan Löfberg
% $Id: callfmincongp.m,v 1.15 2006/07/28 06:27:01 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;

switch options.verbose
    case 0
        options.fmincon.Display = 'off';
    case 1
        options.fmincon.Display = 'final';
    otherwise
        options.fmincon.Display = 'iter';
end

% *********************************
% What type of variables do we have
% *********************************
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Convert to common format for fmincon, mosek and gpposy
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);          

%something failed, perhaps a QP
if problem & isempty(sigmonial_variables) 
    % This is an LP or QP!
    % Go to standard fmincon
    if options.verbose
        disp('Conversion to geometric program failed. Trying general non-convex model in fmincon');
        disp(' ');
    end
    interfacedata.solver.tag = strrep(interfacedata.solver.tag,'-geometric','');
    output = callfmincon(interfacedata);
    return
end

if problem == 0
        
    fmincon_fungp([],[],[],prob); % initialize persistent
    fmincon_congp([],[],[],prob); % initialize persistent

    if isempty(x0)
        x0 = zeros(length(linear_variables),1);
    else
        x0 = x0(linear_variables);
    end

    % Fake logarithm (extend linearly for small values)
    ind = find(x0<1e-2);
    x0(ind) = exp(log(1e-2)+(x0(ind)-1e-2)/1e-2);
    x0 = log(x0);

    % Clean up the bounds (from branch and bound)
    % Note, these bounds are in the
    % logarithmic variables.
    if ~isempty(lb)
        lb = lb(linear_variables);
        ind = find(lb<1e-2);
        lb(ind) = exp(log(1e-2)+(lb(ind)-1e-2)/1e-2);
        lb = log(lb+sqrt(eps));
    end
    if ~isempty(ub)
        ub = ub(linear_variables);
        ind = find(ub<1e-2);
        ub(ind) = exp(log(1e-2)+(ub(ind)-1e-2)/1e-2);
        ub = log(ub+sqrt(eps));
    end

    if options.savedebug
        ops = options.fmincon;
        save fmincongpdebug prob x0 ops lb ub
    end

    solvertime = clock;
    options.fmincon.GradObj    = 'on';
    options.fmincon.GradConstr = 'off';        
    if isfield(options.fmincon,'tlprob')
        [xout,fmin,flag,output,lambda] = fmincon('fmincon_fungp',x0,[],[],[],[],lb,ub,'fmincon_congp',options.fmincon,options.fmincon.tlprob);
    else
        [xout,fmin,flag,output,lambda] = fmincon('fmincon_fungp',x0,[],[],[],[],lb,ub,'fmincon_congp',options.fmincon);
    end
    solvertime = etime(clock,solvertime);
    
    x = zeros(length(c),1);
    x(linear_variables) = exp(xout);
    problem = 0;
    % Check, currently not exhaustive...
    if flag==0
        problem = 3;
    else
        if flag>0
            problem = 0;
        else
            if isempty(x)
                x = repmat(nan,length(c),1);
            end
            if 0%any((A*x-b)>sqrt(eps)) | any( abs(Aeq*x-beq)>sqrt(eps))
                problem = 1; % Likely to be infeasible
            else
                if c'*x<-1e10 % Likely unbounded
                    problem = 2;
                else          % Probably convergence issues
                    problem = 5;
                end
            end
        end
    end
else
    x = [];
    solvertime = 0;
end

% Internal format for duals (currently not supported in GP)
D_struc = [];

infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolverinput
    solverinput.A = [];
    solverinput.b = [];
    solverinput.Aeq = [];
    solverinput.beq = [];
    solverinput.options = options.fmincon;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;