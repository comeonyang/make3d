function output = callquadprog(interfacedata)

% Author Johan Löfberg 
% $Id: callquadprog.m,v 1.12 2006/08/16 15:57:36 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

switch options.verbose
case 0
    options.quadprog.Display = 'off';
case 1
    options.quadprog.Display = 'final';
otherwise
    options.quadprog.Display = 'iter';
end
    
% QUAPROG does not like lb==ub in (LINUX/6.1)
equality_in_bound = find((abs(lb-ub)<1e-12) & ~isinf(lb));
n = length(c);
m = length(equality_in_bound);
if ~isempty(equality_in_bound)
    F_struc = [-lb(equality_in_bound) sparse(1:m,equality_in_bound,ones(m,1),m,n);F_struc];
    ub(equality_in_bound) = ub(equality_in_bound) + 1;
    lb(equality_in_bound) = lb(equality_in_bound) - 1;
    K.f = K.f + m;
end

if ~isempty(F_struc)
    Aeq = -F_struc(1:1:K.f,2:end);
    beq = F_struc(1:1:K.f,1);        
    A =-F_struc(K.f+1:end,2:end);
    b = F_struc(K.f+1:end,1);   
else
    A = [];
    b = [];
    Aeq = [];
    beq = [];
end

if isfield(options.quadprog,'LargeScale')
    if ~isequal(options.quadprog.LargeScale,'on')
        Q = full(Q);
        c = full(c);
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

if options.savedebug
    ops = options.quadprog;
    save quadprogdebug Q c A b Aeq beq lb ub x0 ops
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock; 
if nnz(Q) == 0
    % Quadprog switches anyway, so let us do it to avoid warnings
    [x,fmin,flag,output,lambda] = linprog(c, A, b, Aeq, beq, lb, ub, x0,options.quadprog);
else
    [x,fmin,flag,output,lambda] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub, x0,options.quadprog);
    if flag==5
        [x,fmin,flag,output,lambda] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub, [],options.quadprog);
    end
end
%etime(clock,solvertime)
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

problem = 0;

% Internal format for duals
if ~isempty(lambda)
    D_struc = [lambda.eqlin;lambda.ineqlin];
else
    D_struc = [];
end

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    if flag==-2
        problem = 1;
    else
        if flag>0
            problem = 0;
        else
            if isempty(x)
                x = repmat(nan,length(c),1);
            end
            if any((A*x-b)>sqrt(eps)) | any( abs(Aeq*x-beq)>sqrt(eps))
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
end
infostr = yalmiperror(problem,'QUADPROG');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
    solverinput.beq = beq;
    solverinput.c = c;
    solverinput.H = Q;
    solverinput.options = options.quadprog;
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
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;