function output = callfmincon(interfacedata)

% Author Johan Löfberg
% $Id: callfmincon.m,v 1.33 2006/09/13 12:39:17 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
monomtable = interfacedata.monomtable;

switch options.verbose
    case 0
        options.fmincon.Display = 'off';
    case 1
        options.fmincon.Display = 'final';
    otherwise
        options.fmincon.Display = 'iter';
end

% Do some pre-calc to be used in calls from fmincon
nonlinearindicies = find(interfacedata.variabletype~=0);
nonlinearindicies = union(nonlinearindicies,interfacedata.evalVariables);
linearindicies    = find(interfacedata.variabletype==0);
linearindicies    = setdiff(linearindicies,nonlinearindicies);
interfacedata.nonlinearindicies = nonlinearindicies;
interfacedata.linearindicies    = linearindicies;

% % % Move nonlinear bounds to constraints
% if ~isempty(lb)
%     finite = find(~isinf(lb(nonlinearindicies)));
%     if ~isempty(finite)
%         temp = F_struc(1:K.f,:);
%         F_struc(1:K.f,:) = [];
%         n = length(c);
% 
%         for i = 1:length(finite)
%             j = nonlinearindicies(i);
%             F_struc = [-lb(j) sparse(1,j,1,1,n);F_struc];
%         end
%         K.l = K.l + length(finite);
%         F_struc = [temp;F_struc];
%         interfacedata.K = K;
%         interfacedata.F_struc = F_struc;
%     end
% end
% if ~isempty(ub)
%     finite = find(~isinf(lb(nonlinearindicies)));
%     if ~isempty(finite)
%         temp = F_struc(1:K.f,:);
%         F_struc(1:K.f,:) = [];
%         n = length(c);
% 
%         for i = 1:length(finite)
%             j = nonlinearindicies(i);
%             F_struc = [ub(j) -sparse(1,j,1,1,n);F_struc];
%         end
%         K.l = K.l + length(finite);
%         F_struc = [temp;F_struc];
%         interfacedata.K = K;
%         interfacedata.F_struc = F_struc;
%     end
% end

any_constraints = (K.s+K.f+K.l+K.q)>0;

interfacedata.Anonlinineq = [];
interfacedata.bnonlinineq = [];
interfacedata.Anonlineq = [];
interfacedata.bnonlineq = [];

if K.f>0
    Aeq = -interfacedata.F_struc(1:1:K.f,2:end);
    beq = interfacedata.F_struc(1:1:interfacedata.K.f,1);
        
    nonlinear_equalities_indicies = find(any(Aeq(:,nonlinearindicies),2));
    interfacedata.Anonlineq = Aeq(nonlinear_equalities_indicies,:);
    interfacedata.bnonlineq = beq(nonlinear_equalities_indicies);
       
    Aeq(nonlinear_equalities_indicies,:) = [];
    beq(nonlinear_equalities_indicies,:) = [];
    Aeq(:,nonlinearindicies) = [];
    interfacedata.F_struc(1:interfacedata.K.f,:) = [];
    interfacedata.K.f = 0;   
else
    Aeq = [];
    beq = [];
end

if interfacedata.K.l>0
    A = -interfacedata.F_struc(1:interfacedata.K.l,2:end);
    b = interfacedata.F_struc(1:interfacedata.K.l,1);
         
    nonlinear_inequalities_indicies = find(any(A(:,nonlinearindicies),2));
        
    interfacedata.Anonlinineq = A(nonlinear_inequalities_indicies,:);
    interfacedata.bnonlinineq = b(nonlinear_inequalities_indicies);

    A(nonlinear_inequalities_indicies,:) = [];
    b(nonlinear_inequalities_indicies,:) = [];
    A(:,nonlinearindicies) = [];
    
    interfacedata.F_struc(1:interfacedata.K.l,:) = [];
    interfacedata.K.l = 0;    
else
    A = [];
    b = [];
end

if isfield(options.fmincon,'LargeScale')
    if isequal(options.fmincon.LargeScale,'off')
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

% This helps with robustness in bnb in some cases
x0candidate = zeros(length(c),1);
if ~isempty(lb) & ~isempty(lb)
    bounded = find(~isinf(lb) & ~isinf(ub));
    x0candidate(bounded) = (lb(bounded) + ub(bounded))/2;
    bounded_below = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_below) = lb(bounded_below) + 0.5;
    bounded_above = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_above) = lb(bounded_above) + 0.5;
end
    
if isempty(x0)
    x0 = x0candidate(linearindicies);
else
    if ~isempty(lb) & ~isempty(ub)
        x0((x0 < lb) | (x0 > ub)) = x0candidate((x0 < lb) | (x0 > ub));
    end
    x0 = x0(linearindicies);
end

if ~isempty(lb)
    lb = lb(linearindicies);
end
if ~isempty(ub)
    ub = ub(linearindicies);
end

if options.savedebug
    ops = options.fmincon;
    save fmincondebug interfacedata A b Aeq beq x0 lb ub ops
end

showprogress('Calling FMINCON',options.showprogress);

% Precalc for the callbacks
params = setup_fmincon_params(interfacedata);
if (params.SimpleQuadraticObjective | params.SimpleNonlinearObjective) & isempty(interfacedata.evalMap)
    options.fmincon.GradObj = 'on';
end
params.linearconstraints =  isempty(params.interfacedata.evalMap) & isempty(params.interfacedata.Anonlinineq) & isempty(params.interfacedata.Anonlineq) & isequal( params.interfacedata.K.q,0) & isequal( params.interfacedata.K.s,0);
params.nonlinearinequalities = ~isempty(params.interfacedata.Anonlinineq);
params.nonlinearequalities = ~isempty(params.interfacedata.Anonlineq);

solvertime = clock;
[xout,fmin,flag,output,lambda] = fmincon('fmincon_fun',x0,A,b,Aeq,beq,lb,ub,'fmincon_con',options.fmincon,params);    
solvertime = etime(clock,solvertime);

if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
    x = x(1:length(c));
end

problem = 0;

% Internal format for duals
D_struc = [];

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

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
    solverinput.beq = beq;
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
output = createoutput(x,D_struc,[],problem,'FMINCON',solverinput,solveroutput,solvertime);