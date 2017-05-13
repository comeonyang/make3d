function varargout = entropy(varargin)
%ENTROPY 
%
% y = ENTROPY(x)
%
% Computes/declares entropy -sum(x.*log(x))
%
% Implemented as evalutation based nonlinear operator. Hence, the concavity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

% Author Johan Löfberg
% $Id: entropy.m,v 1.2 2006/08/16 09:06:57 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        varargout{1} = -sum(varargin{1}.*log(varargin{1}));
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects.        
        if min(size(varargin{1}))>1
            error('ENTROPY only defined for vector arguments');
        else
            varargout{1} = yalmip('addEvalVariable',mfilename,varargin{1});
        end
            
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};
                
                % This is different from so called extended operators
                % Just do it!
                F = SetupEvaluationVariable(varargin{:});
                
                % Now add your own code, such as domain constraints
                F = F + set(X > 0);
                
                % Let YALMIP know about convexity etc                
                varargout{1} = F;
                varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','none');
                varargout{3} = X;
                
            case 'milp'
                    varargout{1} = [];
                    varargout{2} = [];
                    varargout{3} = [];                
            otherwise
                error('SDPVAR/LOG called with CHAR argument?');
        end
    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end
