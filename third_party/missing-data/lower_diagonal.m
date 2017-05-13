function M = lower_diagonal(numrows,numcols)

M = ones(numrows,numcols);
for i = 1:numrows
  for j = 1:numcols
    if (i < j) 
      M(i,j) = 0;
    end
  end
end
