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
function [GroundLevel] = CalGroundLevel(defaultPara, ImgInfo, Pair)

% This function generate the ground level at the first pair of process

% =========estimating ground level using ground in imgA and imgB on image A  coordinate

% GroundLevel is in Global_Scale===================
ImgInfo(1).Model.Depth.FitDepth = ImgInfo(1).Model.Depth.FitDepth*Pair.DepthScale(1);
ImgInfo(2).Model.Depth.FitDepth = ImgInfo(2).Model.Depth.FitDepth*Pair.DepthScale(2);
% =================================================
 % cleaning Groundmask for A
RangePercent = 100;
[Dy Dx] = size( ImgInfo(1).Model.Depth.FitDepth);
APositionAll = im_cr2w_cr(ImgInfo(1).Model.Depth.FitDepth, permute(ImgInfo(1).Model.Ray,[2 3 1]));
AY_median = median( APositionAll(2,ImgInfo(1).Model.maskG));
ADistant2Ay_median = (APositionAll(2,ImgInfo(1).Model.maskG) - AY_median);
ANumber_YMedian = round( sum(APositionAll(2,ImgInfo(1).Model.maskG)<AY_median)*RangePercent/100);
[Avalue AIndexSort] = sort(ADistant2Ay_median);
AYmedia_mark = AIndexSort(1:ANumber_YMedian);
AGround_mark = zeros(Dy, Dx);
temp = zeros(sum(ImgInfo(1).Model.maskG(:)) ,1);
temp(AYmedia_mark) = 1;
AGround_mark(ImgInfo(1).Model.maskG) = temp;
AGround_mark = logical(AGround_mark);
   % finishing cleaning Groundmask for A

 % cleaning Groundmask for B
BPositionAll = im_cr2w_cr(ImgInfo(2).Model.Depth.FitDepth, permute(ImgInfo(2).Model.Ray,[2 3 1]));
BWrlPosition = Pair.R'*BPositionAll(:,:)+repmat( -Pair.R'*Pair.T, 1, size(BPositionAll(:,:),2));
BWrlPosition = reshape(BWrlPosition,3,55,[]);
BY_median = median( BWrlPosition(2,ImgInfo(2).Model.maskG));
BDistant2By_median = (BWrlPosition(2,ImgInfo(2).Model.maskG) - BY_median);
BNumber_YMedian = round( sum(BWrlPosition(2,ImgInfo(2).Model.maskG)<BY_median)*RangePercent/100);
[Bvalue BIndexSort] = sort(BDistant2By_median);
BYmedia_mark = BIndexSort(1:BNumber_YMedian);
BGround_mark = zeros(Dy, Dx);
temp = zeros(sum(ImgInfo(2).Model.maskG(:)) ,1);
temp(BYmedia_mark) = 1;
BGround_mark(ImgInfo(2).Model.maskG) = temp;
BGround_mark = logical( BGround_mark);
 % finishing cleaning Groundmask for B

% find the jointly median of the ground of image AB in Y direction
   GroundLevel = median([ APositionAll(2,AGround_mark) BWrlPosition(2,BGround_mark)]); % Now Ground is in Img1 Img2 pair scale
	GroundLevel = GroundLevel/Pair.DepthScale(1); % rescale to the Img1 local scale

return;
