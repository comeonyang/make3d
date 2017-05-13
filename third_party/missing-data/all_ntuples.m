function y = all_ntuples(k,n)
% Pick all k-tuples of numbers less than n.  Order is not important.
% Output is a matrix that is n choose k rows and k columns.
if k > n | k < 2
  if k == 1
    y = (1:n)';
  else 
    y = []
  end
else
  if k == n
    y = 1:n;
  else
    smallertups = all_ntuples(k-1, n-1);
    smallertups_and_n = [(repeat(n,size(smallertups,1)))',smallertups];
    y = [smallertups_and_n', (all_ntuples(k,n-1))']';
  end
end