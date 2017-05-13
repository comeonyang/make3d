function y = manyrands_min(m,k,min_vals)
% Given a matrix, take all k-tuples of the columns.  Then for each row,
% take the min value.  Then take the average of all these.
tot = 0;
n = 0;
for tup = (all_ntuples(k,size(m,2)))'
    if k == 1
      m1 = m(:,tup);
    else
      m1 = (min((m(:,tup'))'))';
    end
    % m1 are values that will be picked.
    tot = tot + sum(m1 == min_vals);
    n = n + size(m,1);
end
y = tot/n;