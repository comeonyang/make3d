function [M1, M2] = split_matrix(M)
% Splits a matrix so the odd and even rows go to two different matrices.
%   This is useful for structure-from-motion.
numrows = size(M,1);
M1 = M(1:2:numrows,:);
M2 = M(2:2:numrows,:);
