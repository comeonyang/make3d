function varargout = sdpfun(varargin)
%SDPFUN Gateway to general functions on SDPVAR variables (overloaded)

% Author Johan Löfberg
% $Id: sdpfun.m,v 1.2 2006/09/08 08:52:12 joloef Exp $

switch class(varargin{1})

    case {'struct','double'} % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        varargout{1} = feval(varargin{end},varargin{1:end-1});

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        % try to figure out size of expected output (many->1 or many->many
        varargin_copy = varargin;
        varargin_copy{1} = zeros(size(varargin_copy{1},1),size(varargin_copy{1},2));
        output = feval(varargin_copy{end},varargin_copy{1:end-1});
        if prod(size(output)) == 1
            % MANY -> 1
            varargout{1} = yalmip('addEvalVariable',mfilename,varargin{:});
        else   
            % MANY -> MANY
            y = [];
            [n,m] = size(varargin{1});
            varargin{1} = varargin{1}(:);
            for i = 1:length(varargin{1})
                varargin_copy = varargin;
                varargin_copy{1} = varargin_copy{1}(i);
                y = [y;yalmip('addEvalVariable',mfilename,varargin_copy{:})];
            end
            varargout{1} = reshape(y,n,m);
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}

            case {'graph','milp'}
                t = varargin{2};
                X = varargin{3};

                % This is different from so called extended operators
                % Just do it!
                F = SetupEvaluationVariable(varargin{:});

                % Let YALMIP know about convexity etc
                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp');
                varargout{3} = X;
            otherwise
                error('SDPVAR/SDPFUN called with CHAR argument?');
        end
    otherwise
        error('SDPVAR/SDPFUN called with CHAR argument?');
end
