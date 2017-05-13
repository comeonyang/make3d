function y = compare_real_motion(fo,D,s,n)
% This is meant to be like compare_motion, but with a real motion matrix
% instead of generating a random one.
num_random_starts = 5;
nframes = size(D,1)/2;
npoints = size(D,2);

for j = 1:n
  INC = motion_incidence(fo,nframes,npoints);
  Dtj = remove_translations(D, INC);
  M = approx_full_matrix(Dtj,3);
  ERR = Dtj - M;
  if rem(j,5)==0
    y = [y', compare(M,INC,3,ERR,s,num_random_starts)']';
  else
    y = [y', compare(M,INC,3,ERR,s,num_random_starts)']';
  end
end
  
% First of all, we remove translation as best we can from the occluded data,
% and then give all the methods the same translation-free problem.  This
% should provide a level playing field, while indicating the noise that
% occurs when there might be errors in determining translation.  The other
% alternative would be to give ground truth a translation based on all the
% data.

% Now, M is our pseudo-ground truth; ie. the best guess to the error free matrix.
% M + ERR is Dtj, which is the data matrix with translation removed, as best
%   we can.  The "error" in ground truth will be the difference between our
%   best guess of M and the data; sounds ok.
% If we feed M and ERR to 'compare', we'll be ok since M is only used:
% a) added to ERR to give the data matrix.  This is in fact what we're trying
% to approximate with a rank three matrix.  and b) as the ground truth.
