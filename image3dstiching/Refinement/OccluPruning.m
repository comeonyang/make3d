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
function [Region1Thre PointPix1Thre POriReprojM1Thre PoccluM1Thre  ...
		Region2Thre PointPix2Thre POriReprojM2Thre PoccluM2Thre ] = ...
                OccluPruning(Mask1, Mask2, Region1, PointPix1, POriReprojM1, PoccluM1, OccluDist1, ...
                             Region2, PointPix2, POriReprojM2, PoccluM2, OccluDist2)

% This function Prune the Occlusion point by how far away the occluded plane is
NumBadDepth1 = sum(~Mask1);
NumBadDepth2 = sum(~Mask2);
Pick1 = randperm(sum(Mask1));
Pick2 = randperm(sum(Mask2));

Region1Thre = Region1;
Region1Thre(:,~Mask1) = zeros(4,NumBadDepth1);
PointPix1Thre = PointPix1;
PointPix1Thre(:,~Mask1) = zeros(2,NumBadDepth1);
POriReprojM1Thre = POriReprojM1;
POriReprojM1Thre(:,~Mask1) = zeros(2,NumBadDepth1);
PoccluM1Thre = PoccluM1;
PoccluM1Thre(:,~Mask1) = zeros(2, NumBadDepth1);

Region2Thre = Region2;
Region2Thre(:,~Mask2) = zeros(4,NumBadDepth2);
PointPix2Thre = PointPix2;
PointPix2Thre(:,~Mask2) = zeros(2,NumBadDepth2);
POriReprojM2Thre = POriReprojM2;
POriReprojM2Thre(:,~Mask2) = zeros(2,NumBadDepth2);
PoccluM2Thre = PoccluM2;
PoccluM2Thre(:,~Mask2) = zeros(2, NumBadDepth2);

return;
