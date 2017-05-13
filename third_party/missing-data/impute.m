function A = impute(M,INC)
% This will implement the imputation method of Tomasi and Kanade.
% We'll use a very simple, restricted plan to choose where to 
% start and what direction to grow in.  We begin with a rectangle
% in the upper left corner that is fully occupied, and grow one 
% row or column at a time.  We assume there is an obvious maximal
% rectangle to begin with, and that we can only grow in one direction
% at each step.  This will be useful for experiments with simple occlusion patterns.
[numrows,numcols] = size(M);

curcol = 0;
i = 1;
while (i <= numcols & INC(1,i)),
  i = i + 1;
end
curcol = i - 1;

currow = 0;
i = 1;
while (i <= numrows & all(INC(i,1:curcol))),
  i = i+1;  
end
currow = i - 1;

curM = M(1:currow,1:curcol);
while (currow < numrows | curcol < numcols)
  extend_right = 0;
  extend_down = 0;
  if (curcol < numcols & sum(INC(1:currow,min(numcols,curcol+1))) > 3)
    extend_right = 1;
  end
  if (currow < numrows & sum(INC(min(numrows,currow+1),1:curcol)) > 3)
    extend_down = 1;
  end
  if ~(extend_down | extend_right)
    error('No extension available');
  end
  [T, P] = tk_affine(curM);
  % Note that if extension can go in either direction, we extend right.
  if extend_right
    v = find(INC(1:currow,curcol+1));
    curM = [curM, T*(T(v,:)\M(v,curcol+1))];
    curcol = curcol+1;
  else
    v = find(INC(currow+1,1:curcol));
    curM = [curM; (M(currow+1,v)/P(:,v))*P];
    currow = currow+1;
  end
end
A = curM;
