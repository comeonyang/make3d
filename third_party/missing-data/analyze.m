function [y,z,w] = analyze(a,fo)
% a is a matrix with experimental results.  Each row indicates the results for
% one matrix.  The first five columns are the iterative algorithm with random
% starting points.  The sixth column is my algorithm.  The seventh, mine
% followed by the iterative, the 8th is ground truth, the ninth is ground truth
% followed by the iterative one.
% fo is the fraction of elements that are expected to be missing.  
%
% y is (weighted) average, z is percentage achieving min_val,
% w is standard dev.
%
a = round(10000.*a)./10000;
ab = abs(a);
fp = 1 - fo;
numexps = size(ab,1);
rank3 = sum(ab(:,6))/(numexps *fp);
rank3_it = sum(ab(:,7))/(numexps *fp);
gt = sum(ab(:,8))/(numexps *fp);
gt_it = sum(ab(:,9))/(numexps *fp);

randmat = ab(:,1:5);
rand1 = manyrands(randmat, 1)/fp;
rand3 = manyrands(randmat, 3)/fp;
rand5 = manyrands(randmat, 5)/fp;

y = [rand1, rand3, rand5, rank3, rank3_it, gt, gt_it];
w = [manyrands_dev(randmat,1)/fp, ...
     manyrands_dev(randmat,3)/fp, ...
     manyrands_dev(randmat,5)/fp, ...
     std(ab(:,6))/fp, ...
     std(ab(:,7))/fp, ...
     std(ab(:,8))/fp, ...
     std(ab(:,9))/fp];
%     manyrands_dev(ab(:,6),1)/fp, ...
 %    manyrands_dev(ab(:,7),1)/fp, ...
  %   manyrands_dev(ab(:,8),1)/fp, ...
   %  manyrands_dev(ab(:,9),1)/fp];

min_vals = (min(ab'))';
z = [manyrands_min(randmat,1,min_vals), manyrands_min(randmat,3,min_vals), ...
      manyrands_min(randmat,5,min_vals), manyrands_min(ab(:,7),1,min_vals), ...
      manyrands_min(ab(:,9),1,min_vals)];

