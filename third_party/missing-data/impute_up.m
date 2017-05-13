function A = impute_up(M,INC)
% This is like impute, except it starts with a rectangle in the
% lower right, and grows up and to the left.
M2 = rot90(rot90(M));
INC2 = rot90(rot90(INC));
A2 = impute(M2,INC2);
A = rot90(rot90(A2));
