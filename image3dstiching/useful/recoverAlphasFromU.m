% *  This code was used in the following articles:
% *  [1] Learning 3-D Scene Structure from a Single Still Image, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on 3D Representation for Recognition (3dRR-07), 2007.
% *      (best paper)
% *  [2] 3-D Reconstruction from Sparse Views using Monocular Vision, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on Virtual Representations and Modeling 
% *      of Large-scale environments (VRML), 2007. 
% *  [3] 3-D Depth Reconstruction from a Single Still Image, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      International Journal of Computer Vision (IJCV), Aug 2007. 
% *  [6] Learning Depth from Single Monocular Images, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      In Neural Information Processing Systems (NIPS) 18, 2005.
% *
% *  These articles are available at:
% *  http://make3d.stanford.edu/publications
% * 
% *  We request that you cite the papers [1], [3] and [6] in any of
% *  your reports that uses this code. 
% *  Further, if you use the code in image3dstiching/ (multiple image version),
% *  then please cite [2].
% *  
% *  If you use the code in third_party/, then PLEASE CITE and follow the
% *  LICENSE OF THE CORRESPONDING THIRD PARTY CODE.
% *
% *  Finally, this code is for non-commercial use only.  For further 
% *  information and to obtain a copy of the license, see 
% *
% *  http://make3d.stanford.edu/publications/code
% *
% *  Also, the software distributed under the License is distributed on an 
% * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
% *  express or implied.   See the License for the specific language governing 
% *  permissions and limitations under the License.
% *
% */
function [alphas] = recoverAlphasFromU(R)
%[alphaOut, alphaUp, alphaRight, q] = recoverAlphas(cameraLoc, cameraAt, cameraUp, gripperOut, gripperUp)

%up, right, out

if size(R,2) == 9	%rotation matrix
	
	R = reshape(R, 3,3);
	displayFlag = 0;

% Calculates the rotation about the three camera vectors: out, up, and right

R = R';  % Makes rotation from image frame to gripper frame.

q       = [(sqrt( max( 0, 1 + R(1,1) + R(2,2) + R(3,3) ) ) / 2); 
           (sqrt( max( 0, 1 + R(1,1) - R(2,2) - R(3,3) ) ) / 2); 
           (sqrt( max( 0, 1 - R(1,1) + R(2,2) - R(3,3) ) ) / 2); 
           (sqrt( max( 0, 1 - R(1,1) - R(2,2) + R(3,3) ) ) / 2)]; 
signMat = -1 + 2*(0 <= [1                ;
                        (R(3,2) - R(2,3));
		                (R(1,3) - R(3,1));
		                (R(2,1) - R(1,2))]);
q = q .* signMat;

test = q(2)*q(3) + q(4)*q(1);
if (test > 0.499) %% singularity at north pole
	alphaOut = -2 * atan2(q(2),q(1));
	alphaUp = -pi/2;
	alphaRight = 0;
elseif (test < -0.499) %% singularity at south pole
	alphaOut = 2 * atan2(q(2),q(1));
	alphaUp = pi/2;
	alphaRight = 0;
else
	sqx = q(2)*q(2);    sqy = q(3)*q(3);    sqz = q(4)*q(4);
	alphaRight = -atan2(2*q(3)*q(1)-2*q(2)*q(4) , 1 - 2*sqy - 2*sqz);
	alphaUp = -asin(2*test);
	alphaOut = -atan2(2*q(2)*q(1)-2*q(3)*q(4) , 1 - 2*sqx - 2*sqz);
	alphas = [alphaOut, alphaUp, alphaRight];
end

elseif size(R,2) == 3
	alphaOut = acos( R(1));
	alphaRight = asin(R(2) / sqrt(1-R(1)^2));
	alphaUp = acos(R(3) / sqrt(1-R(1)^2));

elseif size(R,2) == 1
	alphaOut = zeros(size(R,1),1);
	alphaRight = zeros(size(R,1),1);
	alphaUp = zeros(size(R,1),1);

elseif size(R,2) == 4
	q = R;
	test = q(2)*q(3) + q(4)*q(1);
	if (test > 0.499) %% singularity at north pole
        	alphaOut = -2 * atan2(q(2),q(1));
	        alphaUp = -pi/2;
        	alphaRight = 0;
	elseif (test < -0.499) %% singularity at south pole
        	alphaOut = 2 * atan2(q(2),q(1));
	        alphaUp = pi/2;
        	alphaRight = 0;
	else
        	sqx = q(2)*q(2);    sqy = q(3)*q(3);    sqz = q(4)*q(4);
	        alphaRight = -atan2(2*q(3)*q(1)-2*q(2)*q(4) , 1 - 2*sqy - 2*sqz);
        	alphaUp = -asin(2*test);
	        alphaOut = -atan2(2*q(2)*q(1)-2*q(3)*q(4) , 1 - 2*sqx - 2*sqz);
        	alphas = [alphaOut, alphaUp, alphaRight];
	end

	alphaOut = acos( R(1));
	alphaRight = asin(R(2) / sqrt(1-R(1)^2));
	alphaUp = acos(R(3) / sqrt(1-R(1)^2));

end;

alphas = [alphaOut, alphaUp, alphaRight];
return;


