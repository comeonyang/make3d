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
function [X]=World2LocalCoord(defaultPara, Xw, Geo, Rotation);

% This function convert the position from World Cartesian coordinate
% to Local coordinate

% First step: to Local coordinate which axis align with north/south and east/west
% Notice: (in world Coordinate) Z axis is not directly as a circle to the
% Latitude, so we can not simple correct the Lantitude by rotating the
% angle of the Lantitude

% 1) Longitude correction
Rw2L1 = Build3DRotationM( 0, 0, -Geo.Longitude, {'x','z','y'}, true);
X1 = Rw2L1*Xw;
% 2) Effective Latitude correction
% Theta_y = atan(X1(3)/X1(1));
% Rw2L2 = Build3DRotationM( 0, Theta_y, 0, {'x','z','y'}, false);
Rw2L2 = Build3DRotationM( 0, Geo.Latitude, 0, {'x','z','y'}, true);
X2 = Rw2L2*X1;

% Second step: to specified camera coordinate
% ,which align with camera primal axis
% 5/25 used the full pitch yaw roll info, R = R_roll*R_pitch*R_yaw
% http://www.microstrain.com/pdf/Orientation%20Conversion%20formulas.pdf
RL2Ca = Build3DRotationM( Rotation(2), -Rotation(1), -Rotation(3), {'x','y','z'}, true);
X = RL2Ca*X2;

return;
