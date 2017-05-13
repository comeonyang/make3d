function V = find_rigid(Mapprox,r)
%NOTE:  ASSUMES RANK 3.  WILL NOT WORK FOR OTHER VALUES OF r.

% extract out the tranformation component, in U
% the structure is then in V

[Ui,Si,Vi] = svd(Mapprox);
%V = Vi(:,1:r)*diag(sqrt(diag(Si(1:r,1:r))));
%U = Ui(:,1:r)*diag(sqrt(diag(Si(1:r,1:r))));
V = Vi(:,1:r)*Si(1:r,1:r);
U = Ui(:,1:r);
u_rows = size(U,1);

for i = 1:u_rows
  Y(i,1:3) = U(i,1)*U(i,:);
  Y(i,4:6) = U(i,2)*U(i,:);
  Y(i,7:9) = U(i,3)*U(i,:);
  X(i) = 1;
end

for i = 1 : u_rows/2
  Y(u_rows+i,1:3) = U(2*i-1,1)*U(2*i,:);
  Y(u_rows+i,4:6) = U(2*i-1,2)*U(2*i,:);
  Y(u_rows+i,7:9) = U(2*i-1,3)*U(2*i,:);
  X(u_rows+i) = 0;
end

%constrain A'A to be symmetric
Y = [Y(:,1) Y(:,2)+Y(:,4) Y(:,3)+Y(:,7) Y(:,5) Y(:,6)+Y(:,8) Y(:,9)];

%use least squares to find A'A
X = lscov(Y,X',eye(u_rows+u_rows/2));

% reshape A'A into a 3x3 matrix
Z(1,1:3) = X(1:3)';
Z(2,2:3) = X(4:5)';
Z(1:3,1) = X(1:3);
Z(3,2) = Z(2,3);
Z(3,3) = X(6);

%Z now contains A'A, where A is the affine xform that minimizes the 
%distance to a rigid transform
[U2,S,V2] = svd(Z);

A = U2 * diag(sqrt(diag(S)));
%A now contains the correct affine xform.

%and UA is the closest rigid transform to U
V = V*inv(A)';


