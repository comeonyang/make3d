function R = remove_translations(M, INC)
% I used to do this incorrectly, thinking I could trivially remove translation 
% with missing data.  This function now ignores INC, and only works for cases 
% with no missing data.

INC;
%ignore this var

trans = -(sum(M'))/size(M,2);
transmat = repeat(trans',size(M,2));
R = M + transmat;
