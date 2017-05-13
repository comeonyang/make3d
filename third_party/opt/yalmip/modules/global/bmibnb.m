function output = bmibnb(p)
%BMIBNB          Branch-and-bound scheme for bilinear programs
%
% BMIBNB is never called by the user directly, but is called by YALMIP from
% SOLVESDP, by choosing the solver tag 'bmibnb' in sdpsettings
%
% The behaviour of BMIBNB can be altered using the fields in the field
% 'bmibnb' in SDPSETTINGS
%
% bmibnb.lowersolver    - Solver for lower bound [solver tag ('')]
% bmibnb.uppersolver    - Solver for upper bound [solver tag ('')]
% bmibnb.lpsolver       - Solver for LP bound tightening [solver tag ('')]
% bmibnb.branchmethod   - Branch strategy ['maxvol' | 'best' ('best')]
% bmibnb.branchrule     - Branch position ['omega' | 'bisect' ('omega')]
% bmibnb.lpreduce       - Improve variable bounds using LP [ real [0,1] (0 means no reduction, 1 means all variables)
% bmibnb.lowrank        - Partition variables into two disjoint sets and branch on smallest [ 0|1 (0)]
% bmibnb.target         - Exit if upper bound<target [double (-inf)]
% bmibnb.roottight      - Improve variable bounds in root using full problem [ 0|1 (1)]
% bmibnb.vartol         - Cut tree when x_U-x_L < vartol on all branching variables
% bmibnb.relgaptol      - Tolerance on relative objective error (UPPER-LOWER)/(1+|UPPER|) [real (0.01)]
% bmibnb.absgaptol      - Tolerance on objective error (UPPER-LOWER) [real (0.01)]
% bmibnb.pdtol          - A number is declared non-negative if larger than...[ double (-1e-6)]
% bmibnb.eqtol          - A number is declared zero if abs(x) smaller than...[ double (1e-6)]
% bmibnb.maxiter        - Maximum number nodes [int (100)]
% bmibnb.maxtime        - Maximum CPU time (sec.) [int (3600)]

% Author Johan Löfberg
% $Id: bmibnb.m,v 1.13 2006/09/20 12:43:45 joloef Exp $

bnbsolvertime = clock;
showprogress('Branch and bound started',p.options.showprogress);

% *************************************************************************
% INITIALIZE DIAGNOSTICS ETC
% *************************************************************************
switch max(min(p.options.verbose,3),0)
case 0
    p.options.bmibnb.verbose = 0;
case 1
    p.options.bmibnb.verbose = 1;
    p.options.verbose = 0;
case 2
    p.options.bmibnb.verbose = 2;
    p.options.verbose = 0;
case 3
    p.options.bmibnb.verbose = 2;
    p.options.verbose = 1;
otherwise
    p.options.bmibnb.verbose = 0;
    p.options.verbose = 0;
end

% *************************************************************************
% Assume feasible (this property can be changed in the presolve codes
% *************************************************************************
p.feasible = 1;

% *************************************************************************
% Extract bounds from model
% *************************************************************************
[p.lb,p.ub,used_rows] = findulb(p.F_struc,p.K);

% These two are needed for ex14_1_7
p = fixbounds(p);
p = presolve_bounds_from_equalities(p);
p = preprocess_eval_bounds(p);
p = fixbounds(p);

% *************************************************************************
% The bilinear solver does not support higher order polynomial (i.e. >2). 
% Hence, we need to bilinearize the problem. however, the upper bound
% solver might support polynomial problems, so we should keep a copy of the
% original problem also
% *************************************************************************
n_in = length(p.c);
[p_bilin,changed] = bilinearize(p);
if (p.solver.uppersolver.constraint.equalities.polynomial &  p.solver.uppersolver.objective.polynomial)
    p_bilin.originalModel = p;
    p = p_bilin;
else
    p = p_bilin;
    p.originalModel = p_bilin;
end

% *************************************************************************
% Pre-calc lists of linear/bilinear/nonlinear variables
% *************************************************************************
[p.linears,p.bilinears,p.nonlinears] = compile_nonlinear_table(p);

% *************************************************************************
% Select branch variables. We should branch on all variables that are
% involved in bilinear terms.
% *************************************************************************
p.branch_variables = decide_branch_variables(p);

% *************************************************************************
% Extract bounds from model
% *************************************************************************
p = preprocess_bilinear_bounds(p);
p = preprocess_eval_bounds(p);

% *************************************************************************
% Now reduce the branch variables by removing bilinear terms that only have
% been introduced due to a convex quadratic objective
% *************************************************************************
p = reduce_bilinear_branching_variables(p);

% *************************************************************************
% Is the given initial solution feasible, or is the zero-solution feasible?
% *************************************************************************
[p,x_min,upper] = initializesolution(p);

% *************************************************************************
% Simple pre-solve loop. The loop is needed for variables defined as w =
% x*y, x = t*u,y=..
% *************************************************************************
i = 0;
goon = 1;
while goon
    start = [p.lb;p.ub];
    i = i+1;
    p = updatenonlinearbounds(p);
    p = presolve_bounds_from_equalities(p);
    p = updatenonlinearbounds(p);   
    p = updatenonlinearbounds(p);    
    p = preprocess_eval_bounds(p);
    i = i+1;
    goon = (norm(start-[p.lb;p.ub],'inf') > 1e-1) & i < 50;
    start = [p.lb;p.ub];
end

% *************************************************************************
% Default root bound improvement
% *************************************************************************
%p.lb = (p.lb+1e5)-1e5;
%p.ub = (p.ub+1e5)-1e5;
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
p = root_node_tighten(p,upper);
p = updatenonlinearbounds(p);
p = preprocess_eval_bounds(p);
p = presolve_bounds_from_equalities(p);
p = updatenonlinearbounds(p);

output = yalmip_default_output;
if p.feasible
    
    % *******************************
    % Bounded domain?
    % *******************************
    if ~isempty(p.bilinears)
        involved = unique(p.bilinears(:,2:3));
        if any(isinf(p.lb(involved))) | any(isinf(p.ub(involved)))
            output = yalmip_default_output;
            output.Primal  = x_min;                             
            output.problem = -6;
            output.infostr = yalmiperror(-6);
            output.solved_nodes = 0;           
            return
        end
    end
    
    % *******************************
    % Save time & memory
    % *******************************
    p.options.savesolverinput  = 0;
    p.options.savesolveroutput = 0;
    p.options.saveduals = 0;
    p.options.dimacs = 0;    
    
    % *******************************
    % RUN BILINEAR BRANCH & BOUND
    % *******************************
    [x_min,solved_nodes,lower,upper] = branch_and_bound(p,x_min,upper);
    
    % ********************************
    % CREATE SOLUTION AND DIAGNOSTICS
    % ********************************
    problem = 0;
    if isinf(upper)
        problem = 1;
    end
    if isinf(lower) & (lower<0)
        problem = 2;
    end
    if isinf(lower) & (lower>0)
        problem = 1;
    end    
    if solved_nodes == p.options.bnb.maxiter
        problem = 3;
    end
else
    problem = 1;
    x_min = repmat(nan,length(p.c),1);
    solved_nodes = 0;
end
output = yalmip_default_output;
output.problem = problem;
output.solved_nodes = solved_nodes;
output.Primal       = x_min(1:n_in);
output.infostr      = yalmiperror(output.problem,'BNB');
output.solvertime   = etime(clock,bnbsolvertime);

function model = fixbounds(model);
polynomials = find(model.variabletype > 0);
LU = max([abs(model.lb) abs(model.ub)],[],2);
for i = 1:length(polynomials)
    j = polynomials(i);
    if j<=length(model.lb)
        vars  = find(model.monomtable(j,:));
        bound = 1;
        for k = 1:length(vars)
            bound = bound * LU(vars(k))^model.monomtable(j,vars(k));
        end
        model.lb(j) = max(model.lb(j),-bound);
        model.ub(j) = min(model.ub(j),bound);
    end
end


