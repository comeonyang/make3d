function output = mpcvx(p)
%MPCVX          Approximate multi-parametric programming
%
% MPCVX is never called by the user directly, but is called by
% YALMIP from SOLVESDP, by choosing the solver tag 'mpcvx' in sdpsettings
%
% The behaviour of MPCVX can be altered using the fields
% in the field 'mpcvx' in SDPSETTINGS
%
% WARNING: THIS IS EXPERIMENTAL CODE
%

% Author Johan Löfberg
% $Id: mpcvx.m,v 1.6 2005/05/07 13:53:20 joloef Exp $

% ********************************
% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
mpsolvertime = clock;
showprogress('mpcvx started',p.options.showprogress);

% *******************************
% Display-logics
% *******************************
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

% *******************************
% No reason to save
% *******************************
p.options.saveduals = 0;

% **********************************
% Generate an exploration set
% **********************************
p.solver.subcall = 'callsedumi';
solver    = eval(['@' p.solver.subcall]);    % LP solver
[THETA,problem] = ray_shoot(p,solver);

%plot(polytope(THETA'));hold on;

% **********************************
% Calculate optimal x in each initial vertex
% **********************************
X   = [];
OBJ = [];
for i = 1:size(THETA,2)
    [x_i,obj_i] = solve_node(p,solver,THETA(:,i));
    X = [X x_i];
    OBJ = [OBJ obj_i];
end

% **********************************
% Partition initial set
% **********************************
node_keeper;
node_keeper(THETA,X,OBJ);
T = delaunayn(THETA',{'Qz','Qt'});

% **********************************
% Do algo on all initial simplicies
% **********************************
optimal_simplicies = [];
p.options.mpcvx.eps = 5e-2;
for i = 1:size(T,1)
    optimal_simplicies = [optimal_simplicies mp_simplex(p,solver,T(i,:)')];
end
mpsolvertime = etime(clock,mpsolvertime)

% **********************************
% Create format compatible with MPT
% **********************************
Pn = polytope;
j = 1;
for i = 1:size(optimal_simplicies,2)
    [theta,x,obj] = node_keeper(optimal_simplicies(:,i));
    m = size(theta,1);
    Minv = inv([ones(1,m+1);theta]);

    try
        Pn = [Pn polytope(theta')];
        Gi{j} = x(find(~ismember(1:length(p.c),p.parametric_variables)),:)*Minv(:,1);
        Fi{j} = x(find(~ismember(1:length(p.c),p.parametric_variables)),:)*Minv(:,2:end);        
        j = j + 1;
    catch
        %Gi{j} = polytope([]);
        %Fi{j} = polytope([]);
    end
end

output.problem = problem;
output.Primal            = nan*ones(length(p.c),1);
output.Dual        = [];
output.Slack       = [];
output.infostr      = yalmiperror(output.problem,'MPCVX');
output.solverinput  = 0;
output.solveroutput.Pn = Pn;
output.solveroutput.Fi = Fi;
output.solveroutput.Gi = Gi;
output.solvertime   = mpsolvertime;

function simplex_solution = mp_simplex(p,solver,theta_indicies)

[theta,x,obj] = node_keeper(theta_indicies);
% Parametric dimension
m = size(theta,1);

Minv = inv([ones(1,m+1);theta]);
M1 = Minv(:,1);
M2 = zeros(size(M1,1),length(p.c));
M2(:,p.parametric_variables) = Minv(:,2:end);
p.F_struc = [M1 M2;p.F_struc];
p.K.l = p.K.l + size(M1,1);

Vbar = obj'*Minv;
c = p.c;
c2 = zeros(length(p.c),1);c2(p.parametric_variables) = -Vbar(2:end);
p.c = p.c+c2;p.c;
output = feval(solver,p);
p.c = c;

upper = obj'*Minv*[1;output.Primal(p.parametric_variables)];
lower = p.c'*output.Primal+output.Primal'*p.Q*output.Primal;

eps_CP_S = min(upper-lower,((upper-lower)/(1+lower)));

% Dig deeper?
%plot(polytope(theta'),struct('color',rand(3,1)))
if eps_CP_S > p.options.mpcvx.eps

    thetac = output.Primal(p.parametric_variables);

    [x_i,obj_i] = solve_node(p,solver,thetac);
    new_index = node_keeper(thetac(:),x_i(:),obj_i);

    simplex_solution = [];
    for i = 1:(size(theta,1)+1)
        j = 1:(size(theta,1)+1);
        j(i)=[];
        theta_test = [theta(:,j) thetac];
        if min(svd([ones(1,size(theta_test,2));theta_test]))>1e-4
            simplex_solution =  [simplex_solution mp_simplex(p,solver,[theta_indicies(j);new_index])];
        end
    end

else
    % This simplex constitutes a node, report back
    simplex_solution = theta_indicies(:);
end

function varargout = node_keeper(varargin)

persistent savedTHETA
persistent savedX
persistent savedOBJ

switch nargin
    case 0 % CLEAR
        savedTHETA = [];
        savedX     = [];
        savedOBJ   = [];
    case 3 % Append
        savedTHETA = [savedTHETA varargin{1}];
        savedX     = [savedX varargin{2}];
        savedOBJ   = [savedOBJ varargin{3}];
        varargout{1} = size(savedTHETA,2);
    case 1 % Get data
        varargout{1} = savedTHETA(:,varargin{1});
        varargout{2} = savedX(:,varargin{1});
        varargout{3} = savedOBJ(:,varargin{1});varargout{3} = varargout{3}(:);
    otherwise
        error('!')
end


function [THETA,problem] = ray_shoot(p,solver)
THETA = [];

p_temp = p;
p_temp.c = p_temp.c*0;
p_temp.Q = 0*p_temp.Q;
for i = 1:25
    p_temp.c(p.parametric_variables) = randn(length(p.parametric_variables),1);
    output = feval(solver,p_temp);
    THETA = [THETA output.Primal(p.parametric_variables)];
end
% Select unique and generate center
THETA = unique(fix(THETA'*1e4)/1e4,'rows')';
center = sum(THETA,2)/size(THETA,2);
THETA = [THETA*0.999+repmat(0.001*center,1,size(THETA,2))];
problem = 0;


function [x_i,obj_i] = solve_node(p,solver,theta);
p_temp = p;
p_temp.F_struc(:,1) = p_temp.F_struc(:,1) + p_temp.F_struc(:,1+p.parametric_variables)*theta;
p_temp.F_struc(:,1+p.parametric_variables) = [];

empty_rows = find(~any(p_temp.F_struc(p.K.f+1:p.K.f+p.K.l,2:end),2));
if ~isempty(empty_rows)
    if all(p_temp.F_struc(p.K.f+empty_rows,1)>=-1e-7)
        p_temp.F_struc(p.K.f+empty_rows,:)=[];
        p_temp.K.l = p_temp.K.l - length(empty_rows);
    else
        feasible = 0;
    end
end

x_var = find(~ismember(1:length(p.c),p.parametric_variables));
theta_var = p.parametric_variables;
Q11 = p.Q(x_var,x_var);
Q12 = p.Q(x_var,theta_var);
Q22 = p.Q(theta_var,theta_var);
c1 = p.c(x_var);
c2 = p.c(theta_var);

p_temp.Q = Q11;
p_temp.c = c1+2*Q12*theta;

%p_temp.c(p.parametric_variables) = [];
%p_temp.Q(:,p.parametric_variables) = [];
%p_temp.Q(p.parametric_variables,:) = [];
output = feval(solver,p_temp);

% Recover complete [x theta]
x_i = zeros(length(p.c),1);
x_i(find(~ismember(1:length(p.c),p.parametric_variables))) = output.Primal;
x_i(p.parametric_variables) = theta;
obj_i = x_i'*p.Q*x_i+p.c'*x_i;
