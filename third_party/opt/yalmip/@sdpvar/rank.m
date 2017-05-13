function varargout=rank(varargin)
%RANK (overloaded)

% Author Johan Löfberg
% $Id: rank.m,v 1.6 2006/05/14 11:54:15 joloef Exp $

% *************************************************************************
% This file defines a nonlinear operator for YALMIP
% Rank is a bit non-standard, so don't look here to learn
% % ***********************************************************************
switch class(varargin{1})

    case 'double' 
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/RANK CALLED WITH DOUBLE. Report error')

    case 'sdpvar' 
        varargout{1} = yalmip('addextendedvariable',mfilename,varargin{1});

    case 'char'
        varargout{1} = set([]);
        properties = struct('convexity','none','monotonicity','none','definiteness','none');
        varargout{2} = properties;
        varargout{3} = varargin{3};

    otherwise
        error('Strange type on first argument in SDPVAR/RANK');
end
