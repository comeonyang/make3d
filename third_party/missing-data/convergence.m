% Here is code -- I have added a few comments to help you.
% N is dimension of square matrix, M is rank we want.  alpha is number
% between 0 and 1 that selects entries that are don't care.
% Program generates random [0,1] NxN matrix A1.  Z is generated from another
% similar random matrix whose elements are then set to one where entries 
% exceed alpha, and zero otherwise.  One entries determine the don't cares.
% Hence, about 1-alpha entries are don't care.  A is set to A1 in the care
% states, zero elsewhere.  Initial random entries are chosen for the
% don't care states in the outer for loop in range -50,+50.
% 
% Inner for loop is really to prevent an infinite loop (crummy programming).
% The real test is the break whent he error change is small.  That's also a
% hack for now.
% 
% XX accumulates sets of all final values of don't care values.

N = 5; M = 2;  alpha = 0.9; NUMSTEPS = 1000;
A1 = rand(N,N);
Z = max(0,sign(rand(N,N) - alpha));
A = A1 - A1.*Z;
Zs = reshape(Z,1,N^2);
index = find(Zs == 1)
count = size(index,2);
X = zeros(N,N);
XX = [];
STEPS = [];

for i = 1:20

 Xs = reshape(X,1,N^2);
 Xf = Xs(index);
 Xf = 100*(rand(1,count)-0.5);
 Xs(index) = Xf;
 X = reshape(Xs,N,N);
% X holds random initial starting conditions for unknown values.
 XfA = [];

 for j = 1:NUMSTEPS
  B = A + X;
% B is the matrix
  [U,D,V] = svd(B);
  BA = V(:,1:M)*D(1:M,1:M)*(V(:,1:M))';
% BA is the best fit rank M matrix to B.
  X = BA.*Z;
% X is updated with missing values from this fit.

  Xs = reshape(X,1,N^2);
  Xfold =Xf;
  Xf = Xs(index);
  XfA = [XfA; Xf];
  Error = sqrt(sum((Xf-Xfold).^2));
% Error is change in magnitude of X.

  A - (BA - X)
  if  Error < 1e-10 
    break
  end
 end
 STEPS = [STEPS; j];
 XX = [XX;Xf];
end
% XX
STEPS
XX(find(STEPS ~= NUMSTEPS), :)