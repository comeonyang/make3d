function y = manyrands_dev(m,k)
% Given a matrix, take all k-tuples of the columns.  Then for each row,
% take the min value.  Then take the average of all these.
tot = 0.0;
n = 0;
for tup = (all_ntuples(k,size(m,2)))'
    m1 = m(:,tup');
    if k == 1
       tot = tot + sum(m1);
    else
      tot = tot + sum(min(m1'));
    end
    n = n + size(m,1);
end
ave = tot/n;


tot = 0.0;
n = 0;
for tup = (all_ntuples(k,size(m,2)))'
    m1 = m(:,tup');
    if k == 1
       tot = tot + sum((m1-ave).^2);
    else
      tot = tot + sum((min(m1')-ave).^2);
    end
    n = n + size(m,1);
end
y = sqrt(tot/(n-1));

