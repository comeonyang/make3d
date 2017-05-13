function [R, means, stds] = trans_best(a, stab, fname)
% a contains the results of the new experiments, that contain the new algorithm
% geared for sfm.  

% Each iteration had 12 results.  All but one had no translation in the
% motion, and did not allow for translation in the algorithm.  There
% were: five applications of Shum with random starting points, the old
% method, old + shum, the new method (geared to motion, using all frame
% pairs, but with no translation allowed), new method + shum, new method 
% allowing for translation, ground truth, ground truth + shum.

% I'm going to ignore the old method here, since it seems clearly worse than the new 
% one, and I won't include it in the paper.

R = [];
for marg = [2, 1.1, 1.01];
% We'll look for all results within 10 percent of the best solution.

numrows = size(a,1);

its = abs(a(:,[1:5,9,12]));
minits = (min(its'))';
minits_marg = marg*minits;
itsbest = its < minits_marg*ones(1,7);

ran1 = sum(sum(itsbest(:,1:5)))/(5*numrows);
ran3 = sum(max((itsbest(:,1:3))'))'/numrows;
ran5 = sum(max((itsbest(:,1:5))'))'/numrows;

Rt = [ran1,ran3,ran5,sum(itsbest(:,6:7))/numrows];
R = [R;Rt];
end

% Here I'm going to process results by producing mean and standard deviation.

shum1 = abs(reshape(a(:,1:5),numrows*5,1));
m1 = mean(shum1);
s1 = std(shum1)/sqrt(5*numrows);
shum3 = (min(abs((a(:,1:3))')))';
shum5 = (min(abs((a(:,1:5))')))';
exps = [shum3, shum5, abs(a(:,6:12))];

stabs = [ones(numrows,2), stab(:,6), stab(:,6), stab(:,8), stab(:,8), ...
              stab(:,10), ones(numrows,2)];
% My methods, and iterative refinements of them may be unstable.
stabexps = exps.*stabs
numstab = sum(stabs)
means = [m1, sum(stabexps)./numstab];
stds = [s1, sqrt( sum(stabexps.^2)./numstab - (sum(stabexps)./numstab).^2)./sqrt(numstab)];

%means = [m1,mean(exps)];
%stds = [s1,std(exps)./sqrt(numrows)];

if fname
fid = fopen(fname, 'a');
fprintf(fid, '%% Data written from trans_best.\n\n');

fprintf(fid, 'Mean/(standard deviation) of: 1it, 3it, 5it, old, oldit, new, newit, newtrans, gt, gtit\n\n');
for i = 1:10
  m = means(i);
  s = stds(i);
  fprintf(fid, '%4.8f/(%4.8f)\n', m, s);
end

fprintf(fid, '\n\nFraction times each iterative method converged to right answer.');
fprintf(fid, '\nMethods are 1it, 3it, 5it, newit, gtit');
fprintf(fid, '\nRows show convergence within 100, 10 and 1 percent of best answer.');

fprintf(fid, '\n\n  %1.4f  %1.4f  %1.4f  %1.4f  %1.4f', R(1,1), R(1,2), R(1,3), R(1,4), R(1,5));
fprintf(fid, '\n  %1.4f  %1.4f  %1.4f  %1.4f  %1.4f', R(2,1), R(2,2), R(2,3), R(2,4), R(2,5));
fprintf(fid, '\n  %1.4f  %1.4f  %1.4f  %1.4f  %1.4f', R(3,1), R(3,2), R(3,3), R(3,4), R(3,5));

fclose(fid);
end 