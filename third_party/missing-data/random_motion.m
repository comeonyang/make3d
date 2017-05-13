function S = random_motion(n)
% We will produce an 2nx4 matrix that describes a random motion.  Each row of
% S is (r1 r2 r3 t) where t is a translation and (r1 r2 r3) is a component of
% rotation.  An odd row and the following even row must have orthonormal 
% rotation vectors.
S = [];
for i = 1:n
   S = [S; one_random_rotation];
end
S = [S,rand(2*n,1)];

function S = one_random_rotation()
r1 = random_normal;
r2 = random_normal;
r2p = r2 - sum(r2.*r1)*r1;
% r2p is projected into the plane normal to r1.  It's no longer a unit vector.
r2pn=r2p./(sqrt(sum(r2p.*r2p)));
S = [r1'; r2pn'];