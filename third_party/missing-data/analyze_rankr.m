function [y] = analyze_rankr(a,fo)
% a is a matrix with experimental results.  Each row indicates the results for
% one matrix.  
% fo is the fraction of elements that are expected to be missing.  
%
% y is (weighted) average, z is percentage achieving min_val,
% w is standard dev.
%
a = round(10000.*a)./10000;
ab = abs(a);
fp = 1 - fo;
numexps = size(ab,1);
rank3 = sum(ab(:,1))/(numexps *fp);
rank3_it = sum(ab(:,2))/(numexps *fp);
dwj_rank3 = sum(ab(:,3))/(numexps *fp);
dwj_rank3_it = sum(ab(:,4))/(numexps *fp);

y = [rank3, rank3_it, dwj_rank3, dwj_rank3_it];
