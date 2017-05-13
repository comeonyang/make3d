function err = affine_error(M,pts,r)

  [U, S, V] = svd(M);
  V = V(:,1:r)*S(1:r,1:r);

%  V = find_rigid(Mapp,r);  % returns 3D structure, up to a rigid xform. 

  pts = pts(1:3,:)';
  A = V\pts;

  %may want to take into account the missing points here -- ie, weight
  %points by the number of frames in which they appeared.
 
  err = sum(sum((V*A - pts).^2));

  
