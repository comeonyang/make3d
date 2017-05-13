function [v] = block_invert(A,r,w)

cntr1 = 1;  % assumes no leading rows or cols of zeros in A
            % will need to fix this+++
cntr2 = 1;
while (size(A) > 0)
  [block, A] = get_block(A,r);
  v(cntr1:r+cntr1-1) = block\w(cntr2:cntr2 + size(block,1) - 1);
  cntr1 = cntr1+r;
  cntr2 = cntr2+size(block,1);
end;
v = v';