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
% $Id: gan.m,v 1.1 2006/08/17 07:57:24 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        x = varargin{1};
        a = varargin{2};
        varargout{1} = sum(x(:).*a(:).^(1./x(:)));
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects.    
        X = varargin{1};
        if min(size(X))>1
            error('ENTROPY only defined for vector arguments');
        else
%             y = [];
%             for i = 1:length(X)
%                 y = [y;yalmip('addEvalVariable',mfilename,X(i),varargin{2})];
%             end
%             varargout{1} = y;
                        varargout{1} = yalmip('addEvalVariable',mfilename,varargin{:});
        end
            
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};
                A = varargin{4};
                % This is different from so called extended operators
                % Just do it!
                F = SetupEvaluationVariable(varargin{:});
                
                % Now add your own code, such as domain constraints
                F = F + set(X > 0);
                
                % Let YALMIP know about convexity etc                
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none');
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
