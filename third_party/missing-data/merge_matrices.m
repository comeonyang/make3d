function M = merge_matrices(M1,M2)
% Opposite of split_matrix.  Recombines.
[numrows1, numcols1] = size(M1);
[numrows2, numcols2] = size(M2);
if numrows1 ~= numrows2;
  error('Matrices must have equal number of rows.');
end
if numcols1 ~= numcols2;
  error('Matrices must have equal number of columns.');
end
  
M = ones(numrows1+numrows2,numcols1);
M(1:2:2*numrows1,:) = M1;
M(2:2:2*numrows1,:) = M2;