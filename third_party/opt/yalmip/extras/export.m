function [model,recoverdata,diagnostic,interfacedata] = export(varargin)
%EXPORT  Exports YALMIP problem to solver specific format
%
%   [MODEL,RECOVERYMODEL,DIAGNOSTIC,INTERNAL] = EXPORT(F,h,options) is the common command to
%   export optimization problems of the following kind
%
%    min        h
%    subject to
%            F >(=) 0
%
%
%   The MODEL is exported in the format defined by the solver chosen
%   in the options structure, or automatically chosen by YALMIP.
%
%   If the solver format not is support by EXPORT,the YALMIP model used to
%   call the solver is returned)
%
%   If YALMIP by some reason failed to generate a model, the DIAGNOSTIC 
%   variable will be non-empty.
%
%   The fourth output is the internal model used by YALMIP to communicate
%   with the generalized solver interfaces.
%
%   The RECOVERYMODEL is used to relate a solution of the exported model
%   to the original variables in YALMIP.

% Arrrgh, new format with logdet much better, but we have to
% take care of old code, requires some testing...
varargin = combatible({varargin{:}});
nargin = length(varargin);
% *********************************
% CHECK INPUT
% *********************************
if nargin<1
    help export
    return
else
    F = varargin{1};
    % Check for wrong syntax
    if ~isempty(F) & ~isa(F,'lmi') & ~isa(F,'constraint')
        error('First argument should be a SET object')
    end

    if isa(F,'constraint')
        F = set(F);
    end
end

model = [];
recoverdata = [];
diagnostic = [];
interfacedata = [];

if nargin>=2
    h = varargin{2};
    if isa(h,'double')
        h = [];
    end
    if ~(isempty(h) | isa(h,'sdpvar') | isa(h,'logdet'))
        error('Second argument (the objective function h) should be an sdpvar or logdet object (or empty).');
    end
    if isa(h,'logdet')
        P  = getP(h);
        g  = getgain(h);
        if g>0
            warning('Perhaps you mean -logdet(P)...')
            diagnostic.yalmiptime = etime(clock,yalmiptime);
            diagnostic.solvertime = 0;
            diagnostic.info = yalmiperror(-2,'YALMIP');
            diagnostic.problem = -2;
            return
        end
        h = getcx(h);
        G = lmi(P>0);
        F = F + G;
    else
        P = [];
        G = [];
    end
else
    P = [];
    h = [];
    G = [];
end
    

if nargin>=3
    options = varargin{3};
    if ~(isempty(options) | isa(options,'struct'))
        error('Third argument (options) should be an sdpsettings struct (or empty).');
    end
    if isempty(options)
        options = sdpsettings;
    end
else
    options = sdpsettings;
end
options.solver = lower(options.solver);

if nargin<6
    if isequal(options.solver,'')
        findallsolvers = 1;
    else
        findallsolvers = 0;
    end
else
    findallsolvers = varargin{6};
end

% Just for safety
if isempty(F) & isempty(G)
    F = lmi;
end

% ******************************************
% COMPILE IN GENERALIZED YALMIP FORMAT
% ******************************************
[interfacedata,recoverdata,solver,diagnostic,F] = compileinterfacedata(F,G,P,h,options,findallsolvers);

if ~isempty(diagnostic)   
    model = [];
    recoverdata = [];
    return
end

% Not official yet
if nargin == 5
    model=interfacedata;
    return
end

% ******************************************
% CONVERT
% ******************************************
switch lower(solver.tag)
            
    case {'sedumi-1.05','sedumi-1.1'}
        
        pars = interfacedata.options.sedumi;
        pars.fid = double(interfacedata.options.verbose);
        model.A = -interfacedata.F_struc(:,2:end);
        model.b = -interfacedata.c;
        model.C = interfacedata.F_struc(:,1);
        model.K = interfacedata.K;
        model.pars = pars;
        
    case 'csdp'
        
        pars = interfacedata.options.csdp;
        pars.printlevel=interfacedata.options.verbose;            
        model.At = -interfacedata.F_struc(:,2:end);
        model.b = -interfacedata.c;
        model.C = interfacedata.F_struc(:,1);
        model.K = interfacedata.K;
        model.pars = pars;
        
    case 'dsdp'
        
        [C,A,b,blk] = sedumi2dsdp(interfacedata.F_struc,interfacedata.c,interfacedata.K);
        interfacedata.options.dsdp.dual_quadratic=spalloc(length(interfacedata.c),length(interfacedata.c),0);
        interfacedata.options.dsdp.printyes = (interfacedata.options.verbose>0);
        model.A = A;
        model.C = C;
        model.b = b;
        model.options = interfacedata.options.dsdp

    case 'sdpa'
        
        [mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(interfacedata.F_struc,interfacedata.c,interfacedata.K);
        if interfacedata.options.verbose==0
            interfacedata.options.sdpa.print = 'no';
        else
            interfacedata.options.sdpa.print = 'display';
        end
        [objVal,x,X,Y,INFO]=sdpam(mDIM,nBLOCK,bLOCKsTRUCT,c,F,[],[],[],options.sdpa);
        model.mDIM = mDIM;
        model.nBLOCK = nBLOCK;
        model.bLOCKsTRUCT = bLOCKsTRUCT;
        model.c = c;
        model.F = F;
        model.x0 = [];
        model.X0 = [];
        model.Y0 = [];
        model.OPTIONS = interfacedata.options.sdpa               
        
    case 'sdpt3-3.1'

        % Convert from internal (sedumi-like) format
        [blk,A,C,b,oldKs]=sedumi2sdpt3(interfacedata.F_struc(:,1),-interfacedata.F_struc(:,2:end),-interfacedata.c,interfacedata.K,30);
        interfacedata.options.sdpt3.printyes=double(interfacedata.options.verbose);
        interfacedata.options.sdpt3.expon=interfacedata.options.sdpt3.expon(1);

        model.blk = blk;
        model.A = A;
        model.C = C;
        model.b = b;
        model.ops = interfacedata.options.sdpt3;

    case 'glpk-glpkmex'
        
        n = length(interfacedata.c);
        LB=repmat(-1e6,n,1);   
        UB=repmat(1e6,n,1);         
        SENSE = 1;     
        C = full(interfacedata.c);   
        A =-interfacedata.F_struc(:,2:end);
        B = full(interfacedata.F_struc(:,1));        
        if length(B)==0;
            A = C';
            B = 1e6;
        end
        CTYPE = [repmat('S',interfacedata.K.f,1); repmat('U',interfacedata.K.l,1)];
        VARTYPE = repmat('C',n,1);
        VARTYPE(interfacedata.integer_variables)='I'; 
        interfacedata.options.glpk.msglev = interfacedata.options.verbose;
        if interfacedata.options.glpk.msglev==1
            interfacedata.options.glpk.msglev = 2;
        end
        model.SENSE = SENSE;
        model.C= C;
        model.A = A;
        model.B = B;
        model.CTYPE = CTYPE;
        model.LB = LB;
        model.UB = UB;
        model.VARTYPE = VARTYPE;
        model.PARAM = interfacedata.options.glpk.lpsolver;
        model.SAVE = interfacedata.options.glpk.save;
                
    case 'pensdp-penopt'
        
        penstruct = sedumi2pen(interfacedata.F_struc,interfacedata.K,interfacedata.c,interfacedata.x0);
        ops = struct2cell(interfacedata.options.pensdp);ops = [ops{1:end}];
        penstruct.ioptions = ops(1:8);
        penstruct.foptions = ops(9:end);
        penstruct.ioptions(4) = interfacedata.options.verbose;
        penstruct.ioptions = penstruct.ioptions;
        penstruct.foptions = penstruct.foptions;
        if penstruct.mconstr == 0
            penstruct.msizes = [];
        end
        model.penstruct = penstruct;
        
    case 'penbmi-penopt'
        
        model.penstruct = sedumi2penbmi(interfacedata.F_struc,interfacedata.c,interfacedata.Q,interfacedata.K,interfacedata.monomtable,interfacedata.options,interfacedata.x0);
               
    otherwise
        model = [];
end



function newinputformat = combatible(varargin)

varargin = varargin{1};

classification = 0;
% 0 : Ambigious
% 1 : Old
% 2 : New

% Try some fast methods to determine...
m = length(varargin);
if m==1
    classification = 2;
elseif m>=3 & isstruct(varargin{3})
    classification = 2;
elseif m>=4 & isstruct(varargin{4})
    classification = 1;
elseif m>=2 & isa(varargin{2},'lmi')
    classification = 1;
elseif m>=3 & isa(varargin{3},'sdpvar')
    classification = 1;
elseif m>=2 & isa(varargin{2},'sdpvar') & min(size(varargin{2}))==1
    classification = 2;
elseif m>=2 & isa(varargin{2},'sdpvar') & prod(size(varargin{2}))>=1
    classification = 1;
elseif m>=2 & isa(varargin{2},'logdet')
    classification = 2;
elseif m==2 & isempty(varargin{2})
    classification = 2;
end

if classification==0
    warning('I might have interpreted this problem wrong due to the new input format in version 3. To get rid of this warning, use an options structure');
    classification = 2;
end

if classification==2
    newinputformat = varargin;
else
    newinputformat = varargin;
    P = varargin{2};
    % 99.9% of the cases....
    if isempty(P)
        newinputformat = {newinputformat{[1 3:end]}};
    else
        if isa(P,'lmi')
            P = sdpvar(P);
        end
        if m>=3
            cxP = newinputformat{3}-logdet(P);
            newinputformat{3}=cxP;
        else
            cxP = -logdet(P);
            newinputformat{3}=cxP;
        end
        newinputformat = {newinputformat{[1 3:end]}};
    end
end
