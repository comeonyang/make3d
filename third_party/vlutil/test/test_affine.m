% TEST_AFFINE  Test AFFINE function

X2 = [1;1] ;
X3 = [1;1;1] ;
X4 = [1;1;1;1];

A2 = eye(2) ;
T2 = eye(2,1) ;
A3 = eye(3) ;
T3 = eye(3,1) ;
A4 = eye(4) ;
T4 = eye(4,1) ;

affine(A2,T2,X2) 
affine(A3,T3,X3) 
affine(A4,T4,X4) 

X2 = ones(2) ;
X3 = ones(3) ;
X4 = ones(4) ;

affine(A2,T2,X2,X2) 
affine(A3,T3,X3,X3,X3) 
affine(A4,T4,X4,X4,X4,X4) 
