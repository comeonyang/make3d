function [x_min,solved_nodes,lower,upper] = branch_and_bound(p,x_min,upper)

% *************************************************************************
% Create handles to solvers
% *************************************************************************
lowersolver = p.solver.lowersolver.call; % For relaxed lower bound problem
uppersolver = p.solver.uppersolver.call; % Local nonlinear upper bound
lpsolver    = p.solver.lpsolver.call;    % LP solver for bound propagation

% *************************************************************************
% GLOBAL PROBLEM DATA (these variables are the same in all nodes)
% *************************************************************************
c       = p.c;
Q       = p.Q;
f       = p.f;
K       = p.K;
options = p.options;

% *************************************************************************
% DEFINE UPPER BOUND PROBLEM. Basically just remove the cuts
% *************************************************************************
p_upper = cleanuppermodel(p);

% *************************************************************************
% Active constraints in main model
% 0   : Inactive constraint (i.e. a cut which unused)
% 1   : Active constraint
% inf : Removed constraint  (found to be redundant)
% *************************************************************************
p.InequalityConstraintState = ones(p.K.l,1);
p.InequalityConstraintState(p.KCut.l,1) = 0;
p.EqualityConstraintState = ones(p.K.f,1);

% *************************************************************************
% LPs ARE USED IN  BOX-REDUCTION
% *************************************************************************
p.lpcuts = p.F_struc(1+p.K.f:1:p.K.l+p.K.f,:);
p.cutState = ones(p.K.l,1);
p.cutState(p.KCut.l,1) = 0; % Don't use to begin with

% *************************************************************************
% INITIALITAZION
% *************************************************************************
p.depth = 0;        % depth in search tree
p.dpos  = 0;        % used for debugging
p.lower = NaN;
lower   = NaN;
gap     = inf;
stack   = [];
solved_nodes = 0;
numGlobalSolutions = 0;

% *************************************************************************
% Silly hack to speed up solver calls
% *************************************************************************
p.getsolvertime = 0;

if options.bmibnb.verbose>0
    disp('* Starting YALMIP bilinear branch & bound.');
    disp(['* Upper solver   : ' p.solver.uppersolver.tag]);
    disp(['* Lower solver   : ' p.solver.lowersolver.tag]);
    if p.options.bmibnb.lpreduce
        disp(['* LP solver      : ' p.solver.lpsolver.tag]);
    end
    disp(' Node       Upper      Gap(%)       Lower    Open');
end

t_start = cputime;
go_on  = 1;

while go_on

    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;

    % ********************************************
    % Reduce size of current box (bound tightening)
    % ********************************************
    p = propagatequadratics(p,upper,lower);
    p = domain_reduction(p,upper,lower,lpsolver);
    p = presolve_bounds_from_equalities(p);
    p = preprocess_eval_bounds(p);

    % ********************************************
    % Detect redundant constraints
    % ********************************************
    p = remove_redundant(p);

    % ********************************************
    % SOLVE LOWER AND UPPER
    % ********************************************
    if p.feasible
        [output,cost] = solvelower(p,options,lowersolver);
        
        % Cplex sucks...
        if output.problem == 12            
            pp = p;
            pp.c = pp.c*0;
            [output2,cost2] = solvelower(pp,options,lowersolver);
            if output2.problem == 0
                output.problem = 2;
            else
                output.problem = 1;
            end
        end
                
            

        info_text = '';
        switch output.problem
            case {1,12} % Infeasible
                info_text = 'Infeasible';
                keep_digging = 0;
                cost = inf;
                feasible = 0;

            case 2 % Unbounded (should not happen!)
                cost = -inf;
                x = output.Primal;

            case {0,3,4} % No problems (disregard numerical problems)

                x = output.Primal;

                % UPDATE THE LOWER BOUND
                if isnan(lower)
                    lower = cost;
                end
                if ~isempty(stack)
                    lower = min(cost,min([stack.lower]));
                else
                    lower = min(lower,cost);
                end

                if cost<upper-1e-5

                    z = evaluate_nonlinear(p,x);

                    % Manage cuts etc
                    p = addsdpcut(p,z);
                    p = addlpcuts(p,x);

                    if numGlobalSolutions < p.options.bmibnb.numglobal
                        [upper,x_min,cost,info_text,numGlobalSolutions] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numGlobalSolutions);
                        [upper,x_min,info_text,numGlobalSolutions] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numGlobalSolutions);
                    end
                else
                    keep_digging = 0;
                    info_text = 'Poor bound';
                end
            otherwise
                cost = -inf;
                x = (p.lb+p.ub)/2;
        end
    else
        info_text = 'Infeasible';
        keep_digging = 0;
        cost = inf;
        feasible = 0;
    end
    solved_nodes = solved_nodes+1;

    % ************************************************
    % PRUNE SUBOPTIMAL REGIONS BASED ON UPPER BOUND
    % ************************************************
    if ~isempty(stack)
        [stack,lower] = prune(stack,upper,options,solved_nodes,p);
    end
    if isempty(stack)
        if isinf(cost)
            lower = upper;
        else
            lower = cost;
        end
    else
        lower = min(lower,cost);
    end

    % ************************************************
    % CONTINUE SPLITTING?
    % ************************************************
    if keep_digging & max(p.ub(p.branch_variables)-p.lb(p.branch_variables))>options.bmibnb.vartol
        spliton = branchvariable(p,options,x);
        bounds  = partition(p,options,spliton,x);
        for i = 1:length(bounds)-1
            node = savetonode(p,spliton,bounds(i),bounds(i+1),-1,x,cost,p.EqualityConstraintState,p.InequalityConstraintState,p.cutState);
            node.bilinears = p.bilinears;
            node = updateonenonlinearbound(node,spliton);
            stack = push(stack,node);
        end
        lower = min([stack.lower]);
    end

    % ************************************************
    %  Pick and create a suitable node
    % ************************************************
    [p,stack] = selectbranch(p,options,stack,x_min,upper);

    if isempty(p)
        if ~isinf(upper)
            relgap = 0;
        end
        if isinf(upper) & isinf(lower)
            relgap = inf;
        end
        depth = 0;
    else
        relgap = 100*(upper-lower)/(1+abs(upper));
        depth = p.depth;
    end
    if options.bmibnb.verbose>0
        fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %s  \n',solved_nodes,upper,relgap,lower,length(stack)+length(p),info_text);
    end

    absgap = upper-lower;
    % ************************************************
    % Continue?
    % ************************************************
    time_ok = cputime-t_start < options.bmibnb.maxtime;
    iter_ok = solved_nodes < options.bmibnb.maxiter;
    any_nodes = ~isempty(p);
    relgap_too_big = (isinf(lower) | isnan(relgap) | relgap>100*options.bmibnb.relgaptol);
    absgap_too_big = (isinf(lower) | isnan(absgap) | absgap>options.bmibnb.absgaptol);
    taget_not_met = upper>options.bmibnb.target;
    go_on = taget_not_met & time_ok & any_nodes & iter_ok & relgap_too_big & absgap_too_big ;
end
if options.bmibnb.verbose>0
    if options.bmibnb.verbose;showprogress([num2str(solved_nodes)  ' Finishing.  Cost: ' num2str(upper) ' Gap: ' num2str(relgap) '%'],options.bnb.verbose);end
end

% *************************************************************************
% Stack functionality
% *************************************************************************
function stack = push(stackin,p)
if ~isempty(stackin)
    stack = [p;stackin];
else
    stack(1)=p;
end

function [p,stack] = pull(stack,method,x_min,upper,branch_variables);
if ~isempty(stack)
    switch method
        case 'maxvol'
            for i = 1:length(stack)
                vol(i) = sum(stack(i).ub(branch_variables)-stack(i).lb(branch_variables));
            end
            [i,j] = max(vol);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);

        case 'best'
            [i,j]=min([stack.lower]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);

        otherwise
    end
else
    p =[];
end

function [stack,lower] = prune(stack,upper,options,solved_nodes,p)
if ~isempty(stack)
    toolarge = find([stack.lower]>upper*(1+1e-4));
    if ~isempty(toolarge)
        stack(toolarge)=[];
    end
    if ~isempty(stack)
        indPOS = find(p.c>0);
        indNEG = find(p.c<0);
        LB = [stack.lb];
        UB = [stack.ub];
        LOWER =  p.c([indPOS(:);indNEG(:)])'*[LB(indPOS,:);UB(indNEG,:)];
        toolarge = find(LOWER > upper*(1+1e-4));
        stack(toolarge)=[];
    end
end
if ~isempty(stack)
    lower = min([stack.lower]);
else
    lower = upper;
end

function node = savetonode(p,spliton,bounds1,bounds2,direction,x,cost,EqualityConstraintState,InequalityConstraintState,cutState);
node.lb = p.lb;
node.ub = p.ub;
node.lb(spliton) = bounds1;
node.ub(spliton) = bounds2;
node.lb(p.integer_variables) = ceil(node.lb(p.integer_variables));
node.ub(p.integer_variables) = floor(node.ub(p.integer_variables));
node.lb(p.binary_variables) = ceil(node.lb(p.binary_variables));
node.ub(p.binary_variables) = floor(node.ub(p.binary_variables));

if direction == -1
    node.dpos = p.dpos-1/(2^sqrt(p.depth));
else
    node.dpos = p.dpos+1/(2^sqrt(p.depth));
end
node.depth = p.depth+1;
node.x0 = x;
node.lpcuts = p.lpcuts;
node.lower = cost;
node.InequalityConstraintState = InequalityConstraintState;
node.EqualityConstraintState = EqualityConstraintState;
node.cutState = cutState;

% *************************************
% DERIVE LINEAR CUTS FROM SDPs
% *************************************
function p = addsdpcut(p,x)
if p.K.s > 0
    top = p.K.f+p.K.l+1;
    newcuts = 1;
    newF = [];
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];
        X = reshape(X,n,n);
        [d,v] = eig(X);
        for m = 1:length(v)
            if v(m,m)<0
                for j = 1:length(x)+1;
                    newF(newcuts,j)= d(:,m)'*reshape(p.F_struc(top:top+n^2-1,j),n,n)*d(:,m);
                end
                % max(abs(newF(:,2:end)),[],2)
                newF(newcuts,1)=newF(newcuts,1)+1e-6;
                newcuts = newcuts + 1;
                if size(p.lpcuts,1)>0
                    dist = p.lpcuts*newF(newcuts-1,:)'/(newF(newcuts-1,:)*newF(newcuts-1,:)');
                    if any(abs(dist-1)<1e-3)
                        newF = newF(1:end-1,:);
                        newcuts = newcuts - 1;
                    end
                end
            end
        end
        top = top+n^2;
    end

    if ~isempty(newF)
        % Don't keep all
        m = size(newF,2);
        %  size(p.lpcuts)
        p.lpcuts = [newF;p.lpcuts];
        p.cutState = [ones(size(newF,1),1);p.cutState];
        violations = p.lpcuts*[1;x];
        p.lpcuts = p.lpcuts(violations<0.1,:);

        if size(p.lpcuts,1)>15*m
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            violations = p.lpcuts*[1;x];
            [i,j] = sort(violations);
            %p.lpcuts = p.lpcuts(j(1:15*m),:);
            %p.cutState = lpcuts = p.lpcuts(j(1:15*m),:);
            %p.lpcuts = p.lpcuts(end-15*m+1:end,:);
        end
    end
end

function p = addlpcuts(p,z)
inactiveCuts = find(~p.cutState);
violation = p.lpcuts(inactiveCuts,:)*[1;z];
need_to_add = find(violation < -1e-4);
if ~isempty(need_to_add)
    p.cutState(inactiveCuts(need_to_add)) = 1;
end
inactiveCuts = find(p.InequalityConstraintState == 0 );
violation = p.F_struc(p.K.f+inactiveCuts,:)*[1;z];
need_to_add = find(violation < -1e-4);
if ~isempty(need_to_add)
    p.InequalityConstraintState(inactiveCuts(need_to_add)) = 1;
end


% *************************************************************************
% Strategy for deciding which variable to branch on
% *************************************************************************
function spliton = branchvariable(p,options,x)
% Split if box is to narrow
width = abs(p.ub(p.branch_variables)-p.lb(p.branch_variables));
if (min(width)/max(width) < 0.1) | (size(p.bilinears,1)==0) %
    [i,j] = max(width);%.*p.weight(p.branch_variables));
    spliton = p.branch_variables(j);
else
    res = x(p.bilinears(:,1))-x(p.bilinears(:,2)).*x(p.bilinears(:,3));
    [ii,jj] = sort(abs(res));
    v1 = p.bilinears(jj(end),2);
    v2 = p.bilinears(jj(end),3);

    acc_res1 = sum(abs(res(find((p.bilinears(:,2)==v1) |  p.bilinears(:,3)==v1))));
    acc_res2 = sum(abs(res(find((p.bilinears(:,2)==v2) |  p.bilinears(:,3)==v2))));

    if ~ismember(v2,p.branch_variables) | (acc_res1>acc_res2)
        spliton = v1;
    else
        spliton = v2;
    end
end

% *************************************************************************
% Strategy for diving the search space
% *************************************************************************
function bounds = partition(p,options,spliton,x_min)
x = x_min;
if isinf(p.lb(spliton))
    p.lb(spliton) = -1e6;
end
if isinf(p.ub(spliton))
    p.ub(spliton) = 1e6;
end

switch options.bmibnb.branchrule
    case 'omega'
        if ~isempty(x_min)
            U = p.ub(spliton);
            L = p.lb(spliton);
            x = x(spliton);
            bounds = [p.lb(spliton) 0.5*max(p.lb(spliton),min(x_min(spliton),p.ub(spliton)))+0.5*(p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
        else
            bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
        end
    case 'bisect'
        bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
    otherwise
        bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
end
if isnan(bounds(2)) %FIX
    if isinf(p.lb(spliton))
        p.lb(spliton) = -1e6;
    end
    if isinf(p.ub(spliton))
        p.ub(spliton) = 1e6;
    end
    bounds(2) = (p.lb(spliton)+p.ub(spliton))/2;
end

function [p,stack] = selectbranch(p,options,stack,x_min,upper)
switch options.bmibnb.branchmethod
    case 'maxvol'
        [node,stack] = pull(stack,'maxvol',x_min,upper,p.branch_variables);
    case 'best'
        [node,stack] = pull(stack,'best',x_min,upper);
    otherwise
        [node,stack] = pull(stack,'best',x_min,upper);
end
% Copy node data to p
if isempty(node)
    p = [];
else
    p.depth = node.depth;
    p.dpos = node.dpos;
    p.lb = node.lb;
    p.ub = node.ub;
    p.lower = node.lower;
    p.lpcuts = node.lpcuts;
    p.x0 = node.x0;
    p.InequalityConstraintState = node.InequalityConstraintState;
    p.EqualityConstraintState = node.EqualityConstraintState;
    p.cutState = node.cutState;
end



function [output,cost] = solvelower(p,options,lowersolver)

removeThese = find(p.InequalityConstraintState==inf);
p.F_struc(p.K.f + removeThese,:) = [];
p.K.l = p.K.l - length(removeThese);

removeThese = find(p.EqualityConstraintState==inf);
p.F_struc(removeThese,:) = [];
p.K.f = p.K.f - length(removeThese);

p_with_bilinear_cuts = p;

if ~isempty(p.bilinears)
    p_with_bilinear_cuts.F_struc(1:p.K.f,:)=[];
    p_with_bilinear_cuts = addmcgormick(p_with_bilinear_cuts);
    p_with_bilinear_cuts.F_struc = [p.F_struc(1:p.K.f,:);p_with_bilinear_cuts.F_struc];
end

if ~isempty(p.evalMap)
    p_with_bilinear_cuts = addEvalVariableCuts(p_with_bilinear_cuts);    
end

% **************************************
% SOLVE NODE PROBLEM
% **************************************
if any(p_with_bilinear_cuts.ub+1e-8<p_with_bilinear_cuts.lb)
    output.problem=1;
    cost = inf;
else
    % We are solving relaxed problem (penbmi might be local solver)
    p_with_bilinear_cuts.monomtable = eye(length(p_with_bilinear_cuts.c));

    if p.solver.lowersolver.objective.quadratic.convex
        % Setup quadratic
        for i = 1:size(p.bilinears,1)
            if p_with_bilinear_cuts.c(p.bilinears(i,1))
                p_with_bilinear_cuts.Q(p.bilinears(i,2),p.bilinears(i,3)) = p_with_bilinear_cuts.c(p.bilinears(i,1))/2;
                p_with_bilinear_cuts.Q(p.bilinears(i,3),p.bilinears(i,2)) = p_with_bilinear_cuts.Q(p.bilinears(i,3),p.bilinears(i,2))+p_with_bilinear_cuts.c(p.bilinears(i,1))/2;
                p_with_bilinear_cuts.c(p.bilinears(i,1)) = 0;
            end
        end
        if ~all(eig(full(p_with_bilinear_cuts.Q))>-1e-12)
            p_with_bilinear_cuts.Q = p.Q;
            p_with_bilinear_cuts.c = p.c;
        end
    end

    fixed = p_with_bilinear_cuts.lb >= p_with_bilinear_cuts.ub;
    if isempty(fixed)
        output = feval(lowersolver,p_with_bilinear_cuts);
        cost = output.Primal'*p_with_bilinear_cuts.Q*output.Primal + p_with_bilinear_cuts.c'*output.Primal + p.f;
        % Minor clean-up
        output.Primal(output.Primal<p.lb) = p.lb(output.Primal<p.lb);
        output.Primal(output.Primal>p.ub) = p.ub(output.Primal>p.ub);
    else
        pp = p_with_bilinear_cuts;
        removethese = fixed;
        if ~isempty(p_with_bilinear_cuts.F_struc)
            p_with_bilinear_cuts.F_struc(:,1)=p_with_bilinear_cuts.F_struc(:,1)+p_with_bilinear_cuts.F_struc(:,1+find(fixed))*p_with_bilinear_cuts.lb(fixed);
            p_with_bilinear_cuts.F_struc(:,1+find(fixed))=[];

            rf = find(~any(p_with_bilinear_cuts.F_struc,2));
            rf = rf(rf<=(p_with_bilinear_cuts.K.f + p_with_bilinear_cuts.K.l));
            p_with_bilinear_cuts.F_struc(rf,:) = [];
            p_with_bilinear_cuts.K.f = p_with_bilinear_cuts.K.f - nnz(rf<=p_with_bilinear_cuts.K.f);
            p_with_bilinear_cuts.K.l = p_with_bilinear_cuts.K.l - nnz(rf>p_with_bilinear_cuts.K.f);
        end
        p_with_bilinear_cuts.c(removethese)=[];
        if nnz(p_with_bilinear_cuts.Q)>0
            p_with_bilinear_cuts.c = p_with_bilinear_cuts.c + 2*p_with_bilinear_cuts.Q(find(~removethese),find(removethese))*p_with_bilinear_cuts.lb(removethese);
            p_with_bilinear_cuts.Q(:,find(removethese))=[];
            p_with_bilinear_cuts.Q(find(removethese),:)=[];
        else
            p_with_bilinear_cuts.Q = spalloc(length(p_with_bilinear_cuts.c),length(p_with_bilinear_cuts.c),0);
        end

        if ~isempty(p_with_bilinear_cuts.binary_variables)
            new_bin = [];
            new_var = find(~fixed);
            for i = 1:length(p_with_bilinear_cuts.binary_variables)
                temp = find(p_with_bilinear_cuts.binary_variables(i) == new_var);
                new_bin =  [new_bin temp(:)'];
            end
            p_with_bilinear_cuts.binary_variables = new_bin;
        end
        if ~isempty(p_with_bilinear_cuts.integer_variables)
            new_bin = [];
            new_var = find(~fixed);
            for i = 1:length(p_with_bilinear_cuts.integer_variables)
                temp = find(p_with_bilinear_cuts.integer_variables(i) == new_var);
                new_bin =  [new_bin temp(:)'];
            end
            p_with_bilinear_cuts.integer_variables = new_bin;
        end
        
        p_with_bilinear_cuts.lb(removethese)=[];
        p_with_bilinear_cuts.ub(removethese)=[];
        p_with_bilinear_cuts.x0(removethese)=[];
        p_with_bilinear_cuts.monomtable(:,find(removethese))=[];
        p_with_bilinear_cuts.monomtable(find(removethese),:)=[];        
        output = feval(lowersolver,p_with_bilinear_cuts);
        x=p.c*0;
        x(removethese)=p.lb(removethese);
        x(~removethese)=output.Primal;
        output.Primal = x;
        cost = output.Primal'*pp.Q*output.Primal + pp.c'*output.Primal + p.f;
    end
end

function output = solveupper(p,p_original,x,options,uppersolver)

% The bounds and relaxed solutions have been computed w.r.t to the relaxed
% bilinear model. We only need the original bounds and variables.
p.lb = p.lb(1:length(p_original.c));
p.ub = p.ub(1:length(p_original.c));
x = x(1:length(p_original.c));

p_upper = p_original;

% ...expand the current node just slightly
p_upper.lb = p.lb;
p_upper.ub = p.ub;
fixed = find(abs([p.lb-p.ub]) < 1e-5);
p_upper.lb(fixed) = (p.lb(fixed) +p.ub(fixed) )/2;
p_upper.ub(fixed) = (p.lb(fixed) +p.ub(fixed) )/2;

% Pick an initial point (this can be a bit tricky...)
% Use relaxed point, shifted towards center of box
if all(x<=p.ub) & all(x>=p.lb)
    p_upper.x0 = 0.1*x + 0.9*(p.lb+p.ub)/2;
else
    p_upper.x0 = (p.lb+p.ub)/2;
end
% Shift towards interior for variables with unbounded lower or upper
lbinfbounds = find(isinf(p.lb));
ubinfbounds = find(isinf(p.ub));
p_upper.x0(ubinfbounds) = x(ubinfbounds)+0.01;
p_upper.x0(lbinfbounds) = x(lbinfbounds)-0.01;
ublbinfbounds = find(isinf(p.lb) & isinf(p.ub));
p_upper.x0(ublbinfbounds) = x(ublbinfbounds);

change_these_lb = setdiff(1:length(p.lb),fixed);
change_these_lb = setdiff(change_these_lb,lbinfbounds);
change_these_ub = setdiff(1:length(p.lb),fixed);
change_these_ub = setdiff(change_these_ub,lbinfbounds);

p_upper.lb(change_these_lb) = 0.99*p.lb(change_these_lb)+p_original.lb(change_these_lb)*0.01;
p_upper.ub(change_these_ub) = 0.99*p.ub(change_these_ub)+p_original.ub(change_these_ub)*0.01;
p_upper.lb(isinf(p_original.lb)) = p_upper.lb(isinf(p_original.lb)) - 0.001;
p_upper.ub(isinf(p_original.ub)) = p_upper.ub(isinf(p_original.ub)) + 0.001;

p_upper.options.saveduals = 0;

ub = p_upper.ub ;
lb = p_upper.lb ;

% Remove redundant equality constraints (important for fmincon)
if p_upper.K.f > 0
    Aeq = -p_upper.F_struc(1:1:p_upper.K.f,2:end);
    beq = p_upper.F_struc(1:1:p_upper.K.f,1);
    redundant = find(((Aeq>0).*Aeq*(p_upper.ub-p_upper.lb) - (beq-Aeq*p_upper.lb) <1e-6));   
    p_upper.F_struc(redundant,:) = [];
    p_upper.K.f = p_upper.K.f - length(redundant);
end

% Solve upper bounding problem
p_upper.options.usex0 = 1;

output = feval(uppersolver,p_upper);
% Project into the box (numerical issue)
output.Primal(output.Primal<p_upper.lb) = p_upper.lb(output.Primal<p_upper.lb);
output.Primal(output.Primal>p_upper.ub) = p_upper.ub(output.Primal>p_upper.ub);





% *************************************************************************
% Heuristics from relaxed
% Basically nothing coded yet. Just check feasibility...
% *************************************************************************
function [upper,x_min,cost,info_text,numglobals] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numglobals)
x(p_upper.binary_variables) = round(x(p_upper.binary_variables));
x(p_upper.integer_variables) = round(x(p_upper.integer_variables));

z = evaluate_nonlinear(p_upper,x);

relaxed_residual = constraint_residuals(p_upper,z);

eq_ok = all(relaxed_residual(1:p_upper.K.f)>=-p_upper.options.bmibnb.eqtol);
iq_ok = all(relaxed_residual(1+p_upper.K.f:end)>=p_upper.options.bmibnb.pdtol);

relaxed_feasible = eq_ok & iq_ok;
info_text = '';
if relaxed_feasible
    this_upper = p_upper.f+p_upper.c'*z+z'*p_upper.Q*z;
    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
        x_min = x;
        upper = this_upper;
        info_text = 'Improved solution';
        cost = cost-1e-10; % Otherwise we'll fathome!
        numglobals = numglobals + 1;
    end
end

% *************************************************************************
% Solve local upper bound problem
% *************************************************************************
function [upper,x_min,info_text,numglobals] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numglobals);

output = solveupper(p,p_upper,x,p.options,uppersolver);
output.Primal(p_upper.integer_variables) = round(output.Primal(p_upper.integer_variables));
output.Primal(p_upper.binary_variables) = round(output.Primal(p_upper.binary_variables));

xu = evaluate_nonlinear(p_upper,output.Primal);

upper_residual = constraint_residuals(p_upper,xu);
feasible = ~isempty(xu) & ~any(isnan(xu)) & all(upper_residual(1:p_upper.K.f)>=-p.options.bmibnb.eqtol) & all(upper_residual(1+p_upper.K.f:end)>=p.options.bmibnb.pdtol);

if feasible
    this_upper = p_upper.f+p_upper.c'*xu+xu'*p_upper.Q*xu;
    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
        x_min = xu;
        upper = this_upper;
        info_text = 'Improved solution';
        numglobals = numglobals + 1;
    end
end

% *************************************************************************
% Detect redundant constraints
% *************************************************************************
function p = remove_redundant(p);

b = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
A = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);

redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));

if length(redundant)>1
    p.InequalityConstraintState(redundant) = inf;
end

if p.options.bmibnb.lpreduce
    b = p.lpcuts(:,1);
    A = -p.lpcuts(:,2:end);
    redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));
    if length(redundant)>1
        p.lpcuts(redundant,:) = [];
        p.cutState(redundant) = [];
    end
end

if p.K.f > 0
    b = p.F_struc(1:p.K.f,1);
    A = -p.F_struc(1:p.K.f,2:end);
    s1 = ((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <1e-6);
    s2 = ((-A>0).*(-A)*(p.ub-p.lb) - ((-b)-(-A)*p.lb) <1e-6);
    redundant = find(s1 & s2);
    if length(redundant)>1
        p.EqualityConstraintState(redundant) = inf;
    end
end

% *************************************************************************
% Clean the upper bound model
% Remove cut constraints, and
% possibly unused variables not used
% *************************************************************************
function p = cleanuppermodel(p);

% We might have created a bilinear model from an original polynomial model.
% We should use the original model when we solve the upper bound problem.
p_bilinear = p;
p = p.originalModel;

% Remove cuts
p.F_struc(p.K.f+p.KCut.l,:)=[];
p.K.l = p.K.l - length(p.KCut.l);
n_start = length(p.c);

% Quadratic mode, and quadratic aware solver?
bilinear_variables = find(p.variabletype == 1 | p.variabletype == 2);
if ~isempty(bilinear_variables)
    used_in_c = find(p.c);
    quadraticterms = used_in_c(find(ismember(used_in_c,bilinear_variables)));
    if ~isempty(quadraticterms) & p.solver.uppersolver.objective.quadratic.nonconvex
        usedinquadratic = zeros(1,length(p.c));
        for i = 1:length(quadraticterms)
            Qij = p.c(quadraticterms(i));
            power_index = find(p.monomtable(quadraticterms(i),:));
            if length(power_index) == 1
                p.Q(power_index,power_index) = Qij;
            else
                p.Q(power_index(1),power_index(2)) = Qij/2;
                p.Q(power_index(2),power_index(1)) = Qij/2;
            end
            p.c(quadraticterms(i)) = 0;
        end
    end
end

% Remove SDP cuts
if length(p.KCut.s)>0
    starts = p.K.f+p.K.l + [1 1+cumsum((p.K.s).^2)];
    remove_these = [];
    for i = 1:length(p.KCut.s)
        j = p.KCut.s(i);
        remove_these = [remove_these;(starts(j):starts(j+1)-1)'];
    end
    p.F_struc(remove_these,:)=[];
    p.K.s(p.KCut.s) = [];
end
p.lb = p_bilinear.lb(1:length(p.c));
p.ub = p_bilinear.ub(1:length(p.c));
p.bilinears = [];