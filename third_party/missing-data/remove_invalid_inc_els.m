function [m, inc] = remove_invalid_inc_els(M, INC, r)
% Take subset of M and INC which have all rows and columns with at least r elements.
  valid_cols = find(sum(INC)>r);
  valid_rows = find(sum(INC')>r);
while (size(valid_cols,2) < size(INC,2) | size(valid_rows,2) < size(INC,1))
  M = M(valid_rows,valid_cols);
  INC = INC(valid_rows,valid_cols);
  valid_cols = find(sum(INC)>r);
  valid_rows = find(sum(INC')>r);
end
m=M;
inc=INC;