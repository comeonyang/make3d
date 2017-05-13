function [error, Mapprox, stable] = rankr(M,INC,r,n, iter)
% Fit a rank r matrix to a matrix with missing data.  This uses my
% algorithm of taking the most orthogonal space to the null space generated
% by r-tuples of the columns of M.  In principle, it would be nice
% to use all r-tuples, but instead we'll pick n random ones.


%global NULLSPACEexternal

NULLSPACE = create_nullspace(M,INC,r,n);

if size(NULLSPACE) == [0 0]
  error = -2;
  Mapprox = [];
else
  [error, Mapprox, stable1, stable2] = approximate(M,INC,r,NULLSPACE);
  stable = stable1 & stable2;

% NOTE:  needed to comment out this next re-running step for intensity
% simulations because matrix was too large, so matlab ran out of memory.

%  if (~ stable2) | (error == -2)
%    fprintf(1, 'Repeating procedure....\n');
%    NULLSPACE = create_nullspace(M',INC',r,n);
%    if size(NULLSPACE) == [0 0]
%      error2 = -2;
%    else
%      [error2, Mapprox2, stable1, stable2] = approximate(M',INC',r,NULLSPACE);
%
%     if ((error2 < error) & (error2 ~= -2)) | (error == -2)
%        fprintf(1,'Using new results.\n');
%        error = error2;
%        Mapprox = Mapprox2';
%        stable = stable1 & stable2;
%      end    
%    end
%  end

end
%save data if error.
%if ((error > 3) & stable)
%  fname = sprintf('error%d.bin',iter);
%  fid = fopen(fname, 'W');
%  count = fwrite(fid,[size(M)], 'float');
%  count = fwrite(fid,[M], 'float');
%  count = fwrite(fid,[INC], 'float');
%  count = fwrite(fid,[size(NULLSPACE)], 'float');
%  count = fwrite(fid,[NULLSPACE], 'float');
%  fclose(fid);
%end
%NULLSPACEexternal = NULLSPACE;

if error == -2
  stable = 0;
end
