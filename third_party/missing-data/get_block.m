function [block, B] = get_block(A, r)

[m,n] = size(A);   % number of rows,cols
row = 1;
z = zeros(1,r);
while (A(row,1:r) ~= z)
	row = row + 1;
	if (row > m)
	  break
	end;
end;
block = A(1:row - 1,1:r);
B = A(row:m,r+1:n);
