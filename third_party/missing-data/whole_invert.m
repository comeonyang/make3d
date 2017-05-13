function [U, stable] = whole_invert(V, W, INC, r)
%V is the transformation matrix
%W is the data matrix
%INC is the incidence matrix of the same size as W.  1 denotes valid
%    data and 0 denotes missing data
%r is the desired rank

%to compute each row of U, this function removes the rows of V that
%correspond to missing data, then pre-multiplies the corresponding data
%elements from W by the pseudo-inverse of the new transformation
%matrix. 

i = 1;
n = size(W,1);
W = W';
z = zeros(1,r);
U = zeros(n, r);
stable = 1;

%for each row column of the original W
while (i <= n)

  %remove the rows of V that correspond to missing data
%  sub_mat = V(INC(i,:),:);
% The above line worked in Matlab IV.  The line below replaces it in V.
  sub_mat = V(logical(INC(i,:)),:);

  %check that there is some valid data
%  if (size(sub_mat) > 0)
  if (rank(sub_mat) >= r)

    %choose the valid entries in W and find the entries of U
    %by pre-multiplying by the pseudo-inverse of V
%     t = sub_mat\(W(INC(i,:),i));  % MATLAB IV.
    t = sub_mat\(W(logical(INC(i,:)),i));
  else
    t = z';
    fprintf(1,'WARNING: Rank deficient in whole_invert.\n');
    stable = 0;
%    keyboard
  end;    

  U(i,:) = t';
  i = i + 1;
end;


