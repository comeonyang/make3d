function y = all_rows_sampled(U)
% It's possible that the incidence matrix (and/or bad luck) result
% in our never selecting subsets of the columns of M which have
% entries in a particular row and sufficient rank to provide useful
% information.  If that happens, a column with 1 (or -1?) in that row and
% zero everywhere else will appear in every null space, and, I believe,
% as a component of the svd with a singular value of 0.
last_col = U(:,size(U,2));
y = sum(last_col)~=1;
