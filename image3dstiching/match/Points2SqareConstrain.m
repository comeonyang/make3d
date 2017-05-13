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
function [ Rc, ConS, RoughConS]=Points2SqareConstrain( x, Height)

% This function calculate the Sqare constrain given
% Input:
%	1) Two points
%	2) Heights

% Return:
% Rc2 - (4 x length(x) : the 4 element of each column is the vectorized rotation matrix) 
% ConS - (4 x length(x)) 
%		ConS([1 2],:) - reference corner for the constrain square (x y)
%		ConS([3],:) - sqare width along the epipolar line
%		ConS([4],:) - sqare height othorgonal to the epipolar line
% RoughConS - (4 x length(x))  (not allow rotation)
%		RoughConS([1 2 3 4],:) - [xmax; xmin; ymax; ymin];
xMax = x(1:2,:);
xMin = x(3:4,:);

tRAN = xMax - xMin;
Ptr = tRAN(1,:) < 0;
theta_z = atan(tRAN(2,:)./tRAN(1,:));
theta_z(Ptr) = theta_z(Ptr)+pi;
Rc = [cos(-theta_z); -sin(-theta_z); sin(-theta_z); cos(-theta_z)];

ConS(1:2,:) = xMin;
ConS(3,:) = sum( Rc(1:2,:).*tRAN,1);% also may check sum( Rc(3:4,:).*tRAN,1) close up to zeros to verify correct or not (equal to norm(tRAN,2))
ConS(4,:) = Height;

% generate plot point
tempConS = [ ConS(3,:); ConS(4,:); ...
		zeros(1,size(ConS,2)); ConS(4,:);...
		zeros(1,size(ConS,2)); -ConS(4,:);...
		ConS(3,:); -ConS(4,:)];
Rc_transpose = [cos(theta_z); -sin(theta_z); sin(theta_z); cos(theta_z)];
tempConS(1:2,:) = [sum( Rc_transpose(1:2,:).*tempConS(1:2,:), 1); ...
		      sum( Rc_transpose(3:4,:).*tempConS(1:2,:), 1)];	
tempConS(3:4,:) = [sum( Rc_transpose(1:2,:).*tempConS(3:4,:), 1); ...
		      sum( Rc_transpose(3:4,:).*tempConS(3:4,:), 1)];
tempConS(5:6,:) = [sum( Rc_transpose(1:2,:).*tempConS(5:6,:), 1); ...
		      sum( Rc_transpose(3:4,:).*tempConS(5:6,:), 1)];
tempConS(7:8,:) = [sum( Rc_transpose(1:2,:).*tempConS(7:8,:), 1); ...
		      sum( Rc_transpose(3:4,:).*tempConS(7:8,:), 1)];
tempConS = tempConS + repmat( xMin, 4, 1);

% Rough ConS (not allow rotation) calculated from Plot Points
RoughConS = [ min( tempConS([ 1 3 5 7],:), [], 1); max( tempConS([ 1 3 5 7],:), [], 1);...
                    min( tempConS([ 2 4 6 8],:), [], 1); max( tempConS([ 2 4 6 8],:), [], 1)];


return;
