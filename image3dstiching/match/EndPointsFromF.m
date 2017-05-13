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
function [ EndPointsImg ] = EndPointsFromF(F, x, ScaleImg)

% This function return the two end point on the image for each F*x (epipolar line)

L = F*x; % Epipolar Line (3 by length(x)) each column is a epiplar line vector

% Cut Epipolar Line at image horizontal boundry
x_v_1 = -( 0*L(1,:) + L(3,:))./L(2,:); 		% h_right
x_v_2 = -( ScaleImg(1)*L(1,:) + L(3,:))./L(2,:);  % h_left

% Cut Epipolar Line at image vertical boundry
if x_v_1 < 0
	x_h_1 = -( 0*L(2,:) + L(3,:))./L(1,:); 
	EndPointsImg(1,:) = x_h_1;
	EndPointsImg(2,:) = 0;
elseif x_v_1 > ScaleImg(2)
	x_h_1 = -( ScaleImg(2)*L(2,:) + L(3,:))./L(1,:);
	EndPointsImg(1,:) = x_h_1;
	EndPointsImg(2,:) = ScaleImg(2);
else
	EndPointsImg(1,:) = zeros(1, length(x_v_1));
	EndPointsImg(2,:) = x_v_1;
end 

if x_v_2 < 0
	x_h_2 = -( 0*L(2,:) + L(3,:))./L(1,:); 
	EndPointsImg(3,:) = x_h_2;
	EndPointsImg(4,:) = 0;
elseif x_v_2 > ScaleImg(2)
	x_h_2 = -( ScaleImg(2)*L(2,:) + L(3,:))./L(1,:);
	EndPointsImg(3,:) = x_h_2;
	EndPointsImg(4,:) = ScaleImg(2);
else
	EndPointsImg(3,:) = zeros(1, length(x_v_2));
	EndPointsImg(4,:) = x_v_2;
end 

return;
