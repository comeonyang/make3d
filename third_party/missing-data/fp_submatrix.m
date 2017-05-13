function [submat, subcols] = fp_submatrix(pair,M,INC)
% pair is a 2x1 matrix listing the frames used.  Return the submatrix containing
% the four rows from these frames, and all fully occupied columns, along with
% the indices of these columns.
subrows = [2*pair(1,1)-1, 2*pair(1,1), 2*pair(2,1)-1, 2*pair(2,1)];
subcols = find(sum(INC(subrows,:)) == 4);
submat = M(subrows,subcols);
