function NULLSPACE  = create_nullspace(M,INC,r,n)

% Fit a rank r matrix to a matrix with missing data.  This uses my
% algorithm of taking the most orthogonal space to the null space generated
% by r-tuples of the columns of M.  In principle, it would be nice
% to use all r-tuples, but instead we'll pick n random ones.

%global NULLSPACEexternal

% 10 is the default for the simulated motion experiments.
% however, for the real intensity images, may want to change
% since matrix dimensions are very different, eg, 9 rows
% instead of 40 in M.

%MAXNULLSIZE = 100;  
MAXNULLSIZE = 10;
% if the nullspace gets MAXNULLSIZE times the number of columns in M, we say
% that that's enough, and stop generating elements.

Mdim = size(M);
numrows = Mdim(1);
numcols = Mdim(2);
incsum = sum(INC);
NULLSPACE = [];
current_row = 1;
for i = 1:n
  if size(NULLSPACE,2) >= numrows*MAXNULLSIZE, break, end

  col_nums = find(INC(current_row,:) ~= 0);

% first check that there actually is enough data to recover that row.
% note that although we simply choose columns at random if the check fails,
% it will be impossible to recover the row at any point in the process,
% since it is under-constrained.

  if (length(col_nums) < r)
    col_nums = diff_rand_ints([],r,1,num_cols);
  else
    col_nums = col_nums(diff_rand_ints([],r,1,length(col_nums)));
  end

  incsub = INC(:,col_nums);
  rowsum = sum(incsub')';
  fullrows = rowsum==r;
  if sum(fullrows) > r

    submatrix = M(find(fullrows),col_nums);
    subnull = nulleps(submatrix,.001,-1);
    % .1 is a total hack.  Also, I'm not sure if its good to check sum(fullrows) > r.
    if size(submatrix,1) == size(submatrix,2) + size(subnull,2)
    %% The null space of the submatrix, combined with the submatrix should be a full
    %% rank square matrix.  However, if not, it means that the submatrix was rank
    %% deficient (to within tolerance) and the null matrix is too big.  We can't know
    %% which of these components are the correct ones, and which just arise because the
    %% submatrix was rank deficient.  This test will make the whole method not work if
    %% M really does have rank less than r (in which case the whole null space is right.
      nullTEMP = zeros(numrows,size(subnull,2));
      nullTEMP(fullrows,:) = subnull;
      NULLSPACE = [NULLSPACE, nullTEMP];
    end
  end

  %% Below is a hack added to try to focus attention on the rows for which
  %% we have no information.  Half the time, it makes the current row to be
  %% one for which we don't have any information yet in the nullspace.  The
  %% other half of the time it's just incremented.  I'm not sure whether this
  %% helps significantly.
  if isempty(NULLSPACE)
    sumnonzeros = [];
    nonsamprows = [];
  else
    sumnonzeros = (sum(NULLSPACE' ~= 0)')';
    nonsamprows = find(0 == sumnonzeros);
  end
  if (rand(1) > .5) & (length(nonsamprows) > 0)
    current_row = nonsamprows(random_int(1,length(nonsamprows)));
  else
    current_row = 1+rem(current_row,numrows);
  end
end
%NULLSPACEexternal = NULLSPACE;
