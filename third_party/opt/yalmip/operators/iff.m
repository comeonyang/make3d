function varargout = iff(varargin)
%IFF Logical equivalence
%
% IFF(X,Y) creates a mixed integer representation of
% the constraint X <--> Y, i.e. Y is true iff X is true.
%
% Syntax
%   F = iff(X,Y)
%
% Input
%   X : binary SDPVAR variable or a constraint
%   Y : binary SDPVAR variable or a constraint
%
% Output
%   F : SET object
%
% Examples
%
%  binvar X,Y; F = set(iff(X,Y));
%  sdpvar X;binvar Y; F = set(iff(X>5,Y));
%  sdpvar X;binvar Y; F = set(iff(Y,X==5));
%
% Overloading
%
% The iff overloads == for logic constraints.
%
%  sdpvar X;binvar Y; F = set((X>5) == Y);
%  sdpvar X;binvar Y; F = set(Y == (X==5));

%
% Note
%  The function IFF is not complete, but will be
%  improved upon in future releases.
%
%   See also @SDPVAR/AND, @SDPVAR/OR, IMPLIES

% Author Johan Löfberg
% $Id: iff.m,v 1.2 2006/05/16 07:59:23 joloef Exp $

% There are some cases to take care of...
% X <--> Y  binary/binary
% X <--> Y  binary/(lp,equality)
% X <--> Y  (lp,equality)/binary
% X <--> Y  (lp,equality)/(lp,equality)

X = varargin{1};
Y = varargin{2};

switch class(varargin{1})
    case {'constraint','sdpvar'}
        varargout{1} = set(yalmip('addextendedvariable',mfilename,varargin{:}) == 1);

    case 'char'
        varargout{1} = iff_internal(varargin{3},varargin{4});
        varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp','extra','marker');
        varargout{3} = varargin{3};
end

