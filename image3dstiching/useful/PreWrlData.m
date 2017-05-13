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
function [Position3D ImgPositionWrl] = PreWrlData(K, Depth, I, R, T)

% This function prepare the data for WrlFacestHoriReduce to  use
[VDepth HDepth] = size(Depth);
[VImg HImg] = size(I);
VIndexDepthRes = repmat((1:VDepth)', [1 HDepth]);
HIndexDepthRes = repmat((1:HDepth), [VDepth 1]);
VIndexImgRes = ( VIndexDepthRes -0.5)/VDepth*VImg;
HIndexImgRes = ( HIndexDepthRes -0.5)/HDepth*HImg;
ImgPositionPix = cat(3, HIndexImgRes, VIndexImgRes);
Position3D = inv(K)*reshape( permute(cat(3, ImgPositionPix, ones(VDepth, HDepth)), [3 1 2]), 3, []);
Position3D = Position3D./repmat( sqrt(sum(Position3D.^2, 1)), 3, 1);
Position3D = Position3D.*repmat( Depth(:)', 3, 1);

% rotation and translation
Position3D = R*Position3D+repmat(T, 1, size(Position3D,2));

Position3D = reshape( Position3D, 3, VDepth, []);
Position3D = permute( Position3D, [ 2 3 1]);
Position3D(:,:,3) = -Position3D(:,:,3);

VIndexWrl = (VImg+1-VIndexImgRes - 0.5)/VImg;
HIndexWrl = (HIndexImgRes - 0.5)/HImg;
ImgPositionWrl = cat(3, HIndexWrl, VIndexWrl);
return;

