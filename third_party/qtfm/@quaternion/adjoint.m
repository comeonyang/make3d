function C = adjoint(A, F)
% ADJOINT   Computes the adjoint of the quaternion matrix A.
%
% adjoint(A) or
% adjoint(A, 'complex') returns a complex adjoint matrix.
% adjoint(A, 'real')    returns a real    adjoint matrix.
%
% The definition of the adjoint matrix is not unique (several
% permutations of the layout are possible).  Note that if the
% quaternion A is complexified, it is not possible to compute
% a complex adjoint, since the complex values would have
% complex real and imaginary parts. In this case, the real
% adjoint will work, but its elements will be complex (!).

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout))

if nargin == 1
    F = 'complex'; % Supply the default parameter value.
end

if ~strcmp(F, 'real') & ~strcmp(F, 'complex')
    error('Second parameter value must be ''real'' or ''complex''.')
end

% Extract the components of A. We use scalar() and not s() so that we
% get zero if A is pure.

R = scalar(A); X = x(A); Y = y(A); Z = z(A);

if strcmp(F, 'complex')

    if all(all(imag(R) ~= 0)) | all(all(imag(X) ~= 0)) | ...
       all(all(imag(Y) ~= 0)) | all(all(imag(Z) ~= 0))
        error('Cannot build a complex adjoint with complex elements: use ''real''.');
    end

    % Reference:
    %
    % F. Z. Zhang, Quaternions and Matrices of Quaternions,
    % Linear Algebra and its Applications, 251, January 1997, 21-57.
    
    A1 = complex(R, X); A2 = complex(Y, Z);
    
    C = [[      A1,       A2 ]; ...
         [-conj(A2), conj(A1)]];

else % F must be 'real', since we checked it above.
    
    % Reference:
    %
    % Todd A. Ell, 'Quaternion Notes', 1993-1999, unpublished, defines the
    % layout for a single quaternion. The extension to matrices of quaternions
    % follows easily in similar manner to Zhang above.
    %
    % An equivalent matrix representation for a single quaternion is noted
    % by Ward, J. P., 'Quaternions and Cayley numbers', Kluwer, 1997, p91.
    % It is the transpose of the representation used here.

    C = [[ R,  X,  Y,  Z]; ...
         [-X,  R, -Z,  Y]; ...
         [-Y,  Z,  R, -X]; ...
         [-Z, -Y,  X,  R]];
end

