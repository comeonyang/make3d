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
function [ lamda1 lamda2 Error] = triangulation( defaultPara, R, T, x_calib)

% This function generate the depth from R T and x_calib
if defaultPara.TriLeastSquare
	Q1 = R*x_calib(1:3,:);
	Q2 = x_calib(4:6,:);
	NumRay = size(x_calib,2);
	b = -repmat(T, NumRay, 1);
	A1 = sparse(0,0);
	A2 = sparse(0,0);
	for i=1:NumRay
		A1 = blkdiag(A1, Q1(:,i));	
		A2 = blkdiag(A2, Q2(:,i));	
	end
	lamda = [A1 -A2]\b;
	lamda( lamda <0) = 0;
	lamda1 = lamda(1:NumRay)';
	lamda2 = lamda((NumRay+1):end)';
	Error =  sqrt( sum((R*( x_calib(1:3,:).*repmat(lamda1, 3, 1)) + repmat(T, 1, size(x_calib,2)) - ...
	( x_calib(4:6,:).*repmat(lamda2, 3, 1)) ).^2, 1) );
else
	lamda1 = sdpvar(1,size(x_calib,2));
	lamda2 = sdpvar(1,size(x_calib,2));
	Constrain = set(lamda1 >= 0)+set(lamda2 >= 0);
	sol = solvesdp(Constrain, norm( reshape( R*( x_calib(1:3,:).*repmat(lamda1, 3, 1)) + repmat(T, 1, size(x_calib,2)) - ...
	( x_calib(4:6,:).*repmat(lamda2, 3, 1)), 1, []), 2), defaultPara.opt);
	lamda1 = double(lamda1);
	lamda2 = double(lamda2);

	Error =  sqrt( sum((R*( x_calib(1:3,:).*repmat(lamda1, 3, 1)) + repmat(T, 1, size(x_calib,2)) - ...
	( x_calib(4:6,:).*repmat(lamda2, 3, 1)) ).^2, 1) );
end
return;
