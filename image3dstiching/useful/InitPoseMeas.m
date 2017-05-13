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
function [R T] = InitPoseMeas(defaultPara, ImgInfo1, ImgInfo2);

% This function establish the Pose(Rotation, Translation)
% from the initial measured data (GPS and Compass)
% Input:
% 	defaultPara - camera intrinsic matrix is used
%	ImgInfo - all info of image	
% Return:
% 	R - [R1_2; R2_1] stack rotation matrix of (R1_2) convert local coordinate of image 1 to 2, visa versa
%	T - [T2; T1] stack translation matrix of (T2) translation under local coordinate of image 2, visa versa
 
 % Notice Two Step rotation correctio
 if ~isfield(ImgInfo1,'Rw') || ~isfield(ImgInfo2,'Rw') || ~isfield(ImgInfo1,'Geo') || ~isfield(ImgInfo2,'Geo')
	R = [];
	T = [];
	return; % retrun empty POSE if no Pre Acquire info
 end
 R1_world = Build3DRotationM( ImgInfo1.Rw(2), -ImgInfo1.Rw(1), -ImgInfo1.Rw(3),{'x','y','z'},true)';
 Rworld_2 = Build3DRotationM( ImgInfo2.Rw(2), -ImgInfo2.Rw(1), -ImgInfo2.Rw(3),{'x','y','z'},true);

 % 1) R and T for 1_2
 [X1_2]=World2LocalCoord(defaultPara, ImgInfo1.X_world, ImgInfo2.Geo, ImgInfo2.Rw);
 [X2_2]=World2LocalCoord(defaultPara, ImgInfo2.X_world, ImgInfo2.Geo, ImgInfo2.Rw);
 T2 = (X1_2 - X2_2);
 R1_2 = Rworld_2*R1_world;
%  R1_2_old = Build3DRotationM(-(ImgInfo1.Rw(2)-ImgInfo2.Rw(2)),0,0,{'x','y','z'},true);

 % Important!! permute the entry of T and R according to the image homogeous coordinate
 T2 = T2([2 1 3]);
 R1_2 = R1_2(:,[2 1 3]);
 R1_2 = R1_2([2 1 3],:);
 
 % 2) R and T for 2_1
 [X2_1]=World2LocalCoord(defaultPara, ImgInfo2.X_world, ImgInfo1.Geo, ImgInfo1.Rw);
 [X1_1]=World2LocalCoord(defaultPara, ImgInfo1.X_world, ImgInfo1.Geo, ImgInfo1.Rw);
 T1 = (X2_1- X1_1);
 R2_1 = R1_world'*Rworld_2';
%  R2_1_old = Build3DRotationM(-(ImgInfo2.Rw(2)-ImgInfo1.Rw(2)),0,0,{'x','y','z'},true);

 % Important!! permute the entry of T and R according to the image homogeous coordinate
 T1 = T1([2 1 3]);
 R2_1 = R2_1(:,[2 1 3]);
 R2_1 = R2_1([2 1 3],:);

 % stack
 R = [R1_2; R2_1];
 T = [T2; T1];
 return;
