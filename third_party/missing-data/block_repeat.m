function B = block_repeat(A,p)
% make a sparse, block diagonal matrix by repeating A p times as the blocks.

Adim = size(A);
rows = Adim(1);
cols = Adim(2);

%for i = 1:p
%  block_row = [zeros(cols*(i-1),rows); A'; zeros(cols*(p-i),rows)]';
%  B = [B;block_row];
%end

newrows = rows*p;
newcols = cols*p;

B = sparse([],[],[],newrows,newcols,rows*cols*p);
for i = 0:p-1
  for j = 1:rows
    for k = 1:cols
      B(i*rows+j,i*cols+k) = A(j,k);
    end
  end
end