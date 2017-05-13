function [F,failure,cause] = expandmodel(F,h,options)
% Author Johan Löfberg
% $Id: expandmodel.m,v 1.66 2006/09/13 09:28:51 joloef Exp $

% FIX : Current code experimental, complex, conservative, has issues with
% nonlinearities and is slow...

% All extended variables in the problem. It is expensive to extract this
% one so we will keep it and pass it along in the recursion
extendedvariables = yalmip('extvariables');

% Assume success
failure = 0;
cause = '';

% Early bail out
if isempty(extendedvariables)
    return
end

% Check if it already has ben expanded
already_expanded = expanded(F);
if all(already_expanded)
    if isempty(setdiff(getvariables(h),expanded(h)))
        return
    end
end

% Extract all simple bounds from the model, and update the internal bounds
% in YALMIP. This is done in order to get as tighter big-M models
if ~isempty(F)
    nv = yalmip('nvars');
    yalmip('setbounds',1:nv,repmat(-inf,nv,1),repmat(inf,nv,1));
    LU = getbounds(F);
    yalmip('setbounds',1:nv,LU(:,1),LU(:,2));
end

% Temporary hack to deal with a bug in CPLEX. For the implies operator (and
% some more) YALMIP creates a dummy variable x with set(x==1). Cplex fails
% to solve problem with these stupid variables kept, hence we need to
% remove these variables and constraints...
global MARKER_VARIABLES
MARKER_VARIABLES = [];

% Temporary hack to deal with geometric programs. GPs are messy here,
% becasue we can by mistake claim nonconvexity, since we may have no
% sigmonial terms but indefinite quadratic term, but the whole problem is
% meant to be solved using a GP solver. YES, globals suck, but this is
% only temporary...hrm.
global DUDE_ITS_A_GP
DUDE_ITS_A_GP = 0;

% Keep track of expressions that already have been modelled. Note that if a
% graph-model already has been constructed but we now require a milp, for
% numerical reasons, we should remove the old graph descriptions (important
% for MPT models in particular)
% FIX: Pre-parse the whole problem etc (solves the issues with GP also)
global ALREADY_MODELLED 
global REMOVE_THESE_IN_THE_END 
ALREADY_MODELLED = {};
REMOVE_THESE_IN_THE_END = [];

% Nonlinear operator variables are not allowed to be used in polynomial
% expressions, except if they are exactly modelled, i.e. modelled using
% MILP models. We will expand the model and collect variables that are in
% polynomials, and check in the end if they are exaclty modelled
global OPERATOR_IN_POLYNOM
OPERATOR_IN_POLYNOM = [];

% All variable indicies used in the problem
v1 = getvariables(F);
v2 = depends(F);
v3 = getvariables(h);
v4 = depends(h);

% Speed-hack for LARGE-scale dualizations
if isequal(v3,v4) & isequal(v1,v2)
    variables = uniquestripped([v1 v3]);    
else
    variables = uniquestripped([v1 v2 v3 v4]);
end

% Index to variables modeling operators
extended = find(ismembc(variables,extendedvariables));

if nargin < 3
    options = sdpsettings;
end

% This is a tweak to allow epxansion of bilinear terms in robust problems,
% is expression such as abs(x*w) < 1 for all -1 < w < 1
% This field is set to 1 in robustify and tells YALMIP to skip checking for
% polynomial nonconvexity in the propagation
if ~isfield(options,'expandbilinear')
    options.expandbilinear = 0;
end

% Monomial information. Expensive to retrieve, so we pass this along 
[monomtable,variabletype] = yalmip('monomtable');

% Is this trivially a GP, or meant to be at least?
if strcmpi(options.solver,'gpposy') | strcmpi(options.solver,'fmincon-geometric') | strcmpi(options.solver,'mosek-geometric')
    DUDE_ITS_A_GP = 1;
else
   if ~isequal(options.solver,'fmincon') & ~isequal(options.solver,'') &  ~isequal(options.solver,'mosek')
       % User has specified some other solver, which does not
       % support GPs, hence it cannot be intended to be a GP
       DUDE_ITS_A_GP = 0;
   else
       % Check to see if there are any sigmonial terms on top-level
       DUDE_ITS_A_GP = ~isempty(find(variabletype(variables) == 4));
   end
end

% Constraints generated during recursive process to model operators
F_expand = set([]);

if isempty(F)
    F = set([]);
end

% First, check the objective
variables = uniquestripped([depends(h) getvariables(h)]);
monomtable = monomtable(:,extendedvariables);

% However, some of the variables are already expanded (expand can be called
% sequentially from solvemp and solverobust)
variables = setdiff(variables,expanded(h));

% Determine if we should aim for MILP model directly
if options.allowmilp == 2
    method = 'milp';
else
    method = 'graph';
end

if DUDE_ITS_A_GP == 1
   options.allowmilp = 0;
   method = 'graph'; 
end

% *************************************************************************
% OK, looks good. Apply recursive expansion on the objective
% *************************************************************************
index_in_extended = find(ismembc(variables,extendedvariables));
if ~isempty(index_in_extended)
    extstruct = yalmip('extstruct',variables(index_in_extended));
    if ~isa(extstruct,'cell')
        extstruct = {extstruct};
    end
    [F_expand,failure,cause] = expand(index_in_extended,variables,h,F_expand,extendedvariables,monomtable,'objective',0,options,method,extstruct);
end

% *************************************************************************
% Continue with constraints
% *************************************************************************
constraint = 1;
all_extstruct = yalmip('extstruct');
while constraint <=length(F) & ~failure
    if ~already_expanded(constraint)
        variables = uniquestripped([depends(F(constraint)) getvariables(F(constraint))]);
        [ix,jx,kx] = find(monomtable(variables,:));
        if ~isempty(jx) % Bug in 6.1
            if any(kx>1)
                OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];
            end
        end

        index_in_extended = find(ismembc(variables,extendedvariables));
        if ~isempty(index_in_extended)
            global_index = variables(index_in_extended);
            local_index = [];
            for i = 1:length(global_index)
                local_index = [local_index find(global_index(i) == extendedvariables)];
            end
            extstruct = num2cell(all_extstruct(local_index));
            if is(F(constraint),'equality')
                if options.allowmilp
                    [F_expand,failure,cause] = expand(index_in_extended,variables,-sdpvar(F(constraint)),F_expand,extendedvariables,monomtable,['constraint #' num2str(constraint)],0,options,'milp',extstruct);
                else
                    failure = 1;
                    cause = ['MILP model required for equality in constraint #' num2str(constraint)];
                end
            else
                [F_expand,failure,cause] = expand(index_in_extended,variables,-sdpvar(F(constraint)),F_expand,extendedvariables,monomtable,['constraint #' num2str(constraint)],0,options,method,extstruct);
            end
        end        
    end
    constraint = constraint+1;
end

% *************************************************************************
% Temporary hack to fix the implies operator (cplex has some problem on
% these trivial models where a variable only is used in x==1
% FIX: Automatically support this type of nonlinear operators
% *************************************************************************
if ~isempty(MARKER_VARIABLES)
    MARKER_VARIABLES = sort(MARKER_VARIABLES);
    equalities = find(is(F,'equality'));
    equalities = equalities(:)';
    remove = [];
    for j = equalities
        v = getvariables(F(j));
        if length(v)==1
            if ismembc(v,MARKER_VARIABLES)
                remove = [remove j];
            end
        end
    end
    if ~isempty(remove)
        F(remove) = [];
    end
end

F_expand = lifted(F_expand,1);
% *************************************************************************
% We are done. We might have generated some stuff more than once, but
% luckily we keep track of these mistakes and remove them in the end (this
% happens if we have constraints like set(max(x)<1) + set(max(x)>0) where
% the first constraint would genrate a graph-model but the second set
% creates a milp model.
% *************************************************************************
if ~failure
    F = F + F_expand;    
    if length(REMOVE_THESE_IN_THE_END) > 0
        F = F(find(~ismember(getlmiid(F),REMOVE_THESE_IN_THE_END)));
    end
end

% *************************************************************************
% Normally, operators are not allowed in polynomial expressions. We do
% however allow this if the variable has been modelled with an exact MILP
% model.
% *************************************************************************
Final_model = {ALREADY_MODELLED{unique(OPERATOR_IN_POLYNOM)}};
for i = 1:length(Final_model)
    if ~(strcmp(Final_model{i}.method,'milp') | strcmp(Final_model{i}.method,'none') | options.allownonconvex)
        failure = 1;
        cause = 'Nonlinear operator in polynomial expression.';
        return
    end
end

% declare this model as expanded
F = expanded(F,1);

function [F_expand,failure,cause] = expand(index_in_extended,variables,expression,F_expand,extendedvariables,monomtable,where,level,options,method,extstruct)
global DUDE_ITS_A_GP ALREADY_MODELLED REMOVE_THESE_IN_THE_END OPERATOR_IN_POLYNOM

% *************************************************************************
% Go through all parts of expression to check for convexity/concavity
% First, a small gateway function before calling the recursive stuff
% *************************************************************************
if ~DUDE_ITS_A_GP
    [ix,jx,kx] = find(monomtable(variables,:));
    if ~isempty(jx) % Bug in 6.1
        if any(kx>1)
            OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];          
        end
    end
end

failure = 0;
j = 1;

while j<=length(index_in_extended) & ~failure
    i = index_in_extended(j);
    basis = getbasematrix(expression,variables(i));
    if all(basis >= 0)
        [F_expand,failure,cause] = expandrecursive(recover(variables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,method,extstruct{j},'convex');
    else
        [F_expand,failure,cause] = expandrecursive(recover(variables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,method,extstruct{j},'concave');
    end
    j=j+1;
end