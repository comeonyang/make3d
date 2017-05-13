function nc = percolexp(c)
% given a column of form: [X Y X Y ...]', we generate a set of columns like:
% [X Y Z X Y Z ... ]' which spans the space of all possible 3-D point sets 
% consistent with the images.  These are meant to be many frames, but one pt.
numrows = size(c,1);
numframes = numrows/2;
xs = c(1:2:numrows,1);
ys = c(2:2:numrows,1); 
zs = ones(numframes,1);
nc = ones(1,2*numframes);
for j = 1:2:2*numframes
  nc(j,1) = 
for i = 1:numframes
  
