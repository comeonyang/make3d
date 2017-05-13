function err = trans_affine_error(M,pts)

[U,S,V] = svd(remove_translations(M,0));

V = V(:,1:3)*S(1:3,1:3);
V = [V,ones(size(V,1),1)];

  pts = pts(1:3,:)';
  pts = [pts, zeros(size(pts,1),1)];
  A = V\pts;
  B = V*A;
  err = sum(sum((B(:,1:3) - pts(:,1:3)).^2));
