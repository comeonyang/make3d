function dfdx = jacobian(f,x)
% JACOBIAN Jacobian of scalar or vector
%
% J = JACOBIAN(p)    Jacobian w.r.t all variables in p
% J = JACOBIAN(p,x)  Jacobian w.r.t the SDPVAR variables x
%
% See also SDPVAR, HESSIAN, LINEARIZE

% Author Johan Löfberg
% $Id: jacobian.m,v 1.1 2004/09/15 08:48:21 johanl Exp $

switch nargin    
    case 1        
        dfdx = shadowjacobian(f);
    case 2
        dfdx = shadowjacobian(f,x);
    otherwise
        error('Too many input arguments.');
end
