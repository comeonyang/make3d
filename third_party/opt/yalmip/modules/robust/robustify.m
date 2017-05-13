function [F,h,failure] = robustify(F,h,ops,w)
%ROBUSTIFY  Derives robust counterpart.
%
% [Frobust,objrobust,failure] = ROBUSTIFY(F,h,options) is used to derive
% the robust counterpart of an uncertain YALMIP model.
%
%   min        h(x,w)
%   subject to
%           F(x,w) >(=) 0  for all w in W
%
% The constraints and objective have to satisfy a number of conditions for
% the robustification to be tractable. Please refer to the YALMIP Wiki for
% the current assumptions (this is constantly developing)
%
% See also SOLVEROBUST, UNCERTAIN

% Author Johan Löfberg
% $Id: robustify.m,v 1.19 2006/10/24 12:02:04 joloef Exp $

if nargin < 3
    ops = [];
end

if nargin < 4
    w = [];
end

if isempty(w)
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = recover(getvariables(sdpvar(F(find(unc_declarations)))));
        F = F(find(~unc_declarations));
    else
        error('There is no uncertainty definition in the model.')
    end
end

if isempty(ops)
    ops = sdpsettings;
end

% Figure out which variables are uncertain, certain, and lifted variables
% in the uncertainty description (this code is buggy as ....)
[x,w,x_variables,w_variables,aux_variables,F,failure] = robust_classify_variables(F,h,ops,w);
if failure
    return
end

% Integer variables are OK in x, but not in the uncertainty (robustification
% is based on strong duality in w-space)
integervars = [yalmip('binvariables') yalmip('intvariables')];
ind = find(is(F,'integer') | is(F,'binary'));
if ~isempty(ind)
    integervars = [integervars getvariables(F(ind))];
    if any(ismember(w_variables,integervars))
        failure = 1;
        return
    end
end

% Find  uncertainty description, uncertain and certain constraints
F_w = set([]);
F_x = set([]);
F_xw = set([]);
for i = 1:length(F)
    if all(ismember(depends(F(i)),w_variables))
        % Uncertainty definition
        F_w = F_w + F(i);
    elseif all(ismember(depends(F(i)),x_variables))
        % Certain constraint
        F_x = F_x +  F(i);
    else
        % Uncertain constraint
        F_xw = F_xw + F(i);
    end
end

% Limitation in the modelling language...
if ~isempty(intersect(intersect(depends(F_xw),depends(F_w)),aux_variables))
    disp('You are most likely using a nonlinear operator to describe the');
    disp('uncertainty set (such as norm(w,1) <=1). This is currently not');
    disp('supported. Please model the constraint manually.');
    error('Uncertain model does not satisfy assumptions (nonlinear operator on uncertainty in uncertain constraint)');
end

if length(F_w)==0
    error('There is no uncertainty description in the model.');
end

% Some pre-calc
xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
% Analyze the objective and try to rewrite any uncertainty into the format
% assumed by YALMIP (
if ~isempty(h)
    %
    %[Q,c,f] = quadratic_model(h,xw);
    if 0
        [Q,c,f,dummy,nonquadratic] = quaddecomp(h,xw);
    else
        [Q,c,f,dummy,nonquadratic] = vecquaddecomp(h,xw);
        Q = Q{1};
        c = c{1};
        f = f{1};
    end
    if nonquadratic
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    Q_ww = Q(wind,wind);
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);
 
    if nnz(Q_ww) > 0
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    % Separate certain and uncertain terms, place uncertain terms in the
    % constraints instead
    if is(h,'linear')
        if isempty(intersect(getvariables(w),getvariables(h)))
            h_fixed = h;
        else
            sdpvar t
            F_xw = F_xw + set(h < t);
            h_fixed = t;
            x = [x;t];          
        end
    else
        h_fixed     = x'*Q_xx*x + c_x'*x + f;
        h_uncertain = 2*w'*Q_xw'*x + c_w'*w;
        if ~isa(h_uncertain,'double')
            sdpvar t
            F_xw = F_xw + set(h_uncertain < t);
            h_fixed = h_fixed + t;
            x = [x;t];          
        end
    end
else
    h_fixed = [];
end

% Convert quadratic constraints in uncertainty model to SOCPs.
F_w = convertquadratics(F_w);

% Export uncertainty model to numerical format
ops.solver = '';
[aux1,aux2,aux3,Zmodel] = export(F_w,[],ops,[],[],1);

if ~isempty(Zmodel)
    if length(Zmodel.c) ~= length(w)
        error('Some uncertain variables are unconstrained.')
    end
else
    error('Failed when exporting a model of the uncertainty.')    
end

% OK, we are done with the initial analysis of the involved variables, and
% check of the objective function. 
% 
% At this point, we apply algorithms to robustify constraints (currently we
% only have code for the uncertain conic LP case and polytopic SDP)

F_robust = set([]);

% Pick out the uncertain linear equalities and robustify
F_lp = F_xw(find(is(F_xw,'elementwise')));
F_xw = F_xw - F_lp;
F_robust = F_robust + robustify_lp_conic(F_lp,Zmodel,x,w);

% Pick out uncertain SOCP & SDP constraints and robustify
F_sdp = F_xw(find(is(F_xw,'sdp') | is(F_xw,'socc')));
F_xw = F_xw - F_sdp;
F_robust = F_robust + robustify_sdp_conic(F_sdp,Zmodel,x,w);

% Pick out the uncertain equalities and robustify
F_eq = F_xw(find(is(F_xw,'equality')));
F_xw = F_xw - F_eq;
F_robust = F_robust + robustify_eq_conic(F_eq,Zmodel,x,w);


if length(F_xw) > 0
    error('There are some uncertain constraints that not are supported by YALMIP')
end

% Return the robustfied model
F = F_robust+F_x;
h = h_fixed;

% The model has been expanded, so we have to remember this (trying to
% expand an expanded model leads to nonconvexity error)
F = expanded(F,1); % This is actually done already in expandmodel
h = expanded(h,1); % But this one has to be done manually