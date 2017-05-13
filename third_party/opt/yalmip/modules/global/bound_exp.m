function [Ax, Ay, b] = bound_exp(xL,xU)

% Two lower bounds from tangents
% y > exp(xL) + (x-xL)*exp(xL)
% y > exp(xU) + (x-xU)*exp(xU)

% Upper bound from conneting extreme points
% y < exp(xU)(x-xL)/(xU-xL) +  exp(xL)(xU-x)/(xU-xL)

% can be wrtitten as
% Ax*x + Ay*y < b

Ay = [-1;-1;1];
b = [-exp(xL)+xL*exp(xL);
    -exp(xU)+xU*exp(xU);
    exp(xU)*(-xL)/(xU-xL) +  exp(xL)*(xU)/(xU-xL)];
Ax = [exp(xL);exp(xU); -exp(xU)/(xU-xL) + exp(xL)/(xU-xL)];
