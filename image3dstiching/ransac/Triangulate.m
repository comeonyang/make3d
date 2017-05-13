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
function [lamda1, lamda2] = Triangulate( R, T, X)

% This function triangulate the depths (lamda)
% given camera pose (R T: unit length) of a pair of images 
% Input
%      R - rotation
%      T - unit length translation vector
%      x - calibrated matches point in both images
% Return
%      lamda - triangulated depth for unit length T

X1 = X(1:3,:);
X2 = X(4:6,:);
T = T./(norm(T));
NumMatches = size(X1,2);
M = sparse(0,0);
LastM = [];
X2M = sparse(0,0);
for i = 1:NumMatches
            X2_hat = [[ 0 -X2(3,i) X2(2,i)];...
                      [ X2(3,i) 0 -X2(1,i)];...
              [ -X2(2,i) X2(1,i) 0]];
        LastM = [LastM; X2_hat*T];
	M = blkdiag(M, X2_hat*R*X1(:,i));
	X2M = blkdiag(X2M, X2(:,i));
end
M = [M LastM];
[U S V] =svds(M);
lamda = V(:,end);
lamda = lamda./(lamda(end));
lamda1 = lamda(1:(end-1));
lamda2 = X2M\reshape( R*X1.*repmat(lamda1',3,1)+ repmat(T,1,int32(NumMatches)),[],1);

return;edit Triangulate[
