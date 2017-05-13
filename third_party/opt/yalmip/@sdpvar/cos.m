function varargout = cos (varargin)
%COS (overloaded)

% Author Johan Löfberg
% $Id: cos.m,v 1.1 2006/09/08 08:40:55 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/COS CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if length(varargin{1}) == 1
            varargout{1} = yalmip('addEvalVariable',mfilename,varargin{1});
        else
            y = [];
            for i = 1:length(varargin{1})
                y = [y;yalmip('addEvalVariable',mfilename,extsubsref(varargin{1},i))];
            end
            varargout{1} = y;
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
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;    
                
            otherwise
                error('SDPVAR/COS called with CHAR argument?');
        end
    otherwise
        error('SDPVAR/COS called with CHAR argument?');
end
