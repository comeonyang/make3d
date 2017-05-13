% Test code for the discrete quaternion Fourier transform.
% This code tests the following functions:
%
%  qdft  qdft2
% iqdft iqdft2
%
% It also verifies indirectly many of the basic quaternion operations
% since the qdft depends on them.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Since the two-dimensional qdft code calls the one-dimensional
% code, we need only to test the two dimensionsal code in what
% follows.

% Define one real and one complex quaternion array.

q = quaternion(randn(10,10), randn(10,10), randn(10,10), randn(10,10));
b = quaternion(complex(randn(10,10)),complex(randn(10,10)),...
               complex(randn(10,10)),complex(randn(10,10)));
T = 1e-12;

RA =        unit(quaternion(1,1,1));  % Real axis.
CA = complex(RA, quaternion(1,0,-1)); % Complex axis.

if  isreal(CA) error('Complex axis is not complex.'); end
if ~isreal(RA) error('Real axis is complex.'); end

% Test 1. Verify correct transform and inverse for a real quaternion
% array with a real quaternion axis.

compare(q, iqdft2(qdft2(q, RA, 'L'), RA, 'L'), T, 'qdft failed test 1L.');
compare(q, iqdft2(qdft2(q, RA, 'R'), RA, 'R'), T, 'qdft failed test 1R.');

% Test 2. Verify correct transform and inverse for a real quaternion
% array with a complex axis.

compare(q, iqdft2(qdft2(q, CA, 'L'), CA, 'L'), T, 'qdft failed test 2L.');
compare(q, iqdft2(qdft2(q, CA, 'R'), CA, 'R'), T, 'qdft failed test 2R.');

% Test 3. Verify correct transform and inverse for a complex quaternion
% array with a complex axis.

compare(b, iqdft2(qdft2(b, CA, 'L'), CA, 'L'), T, 'qdft failed test 3L.');
compare(b, iqdft2(qdft2(b, CA, 'R'), CA, 'R'), T, 'qdft failed test 3R.');

% Test 4. Verify correct transform and inverse for a complex quaternion
% array with a real axis.

compare(b, iqdft2(qdft2(b, RA, 'L'), RA, 'L'), T, 'qdft failed test 4L.');
compare(b, iqdft2(qdft2(b, RA, 'R'), RA, 'R'), T, 'qdft failed test 4R.');

clear q b T RA CA
