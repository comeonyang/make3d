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
function [Region1 PointPix1 POriReprojM1 PoccluM1 Region2 PointPix2 POriReprojM2 PoccluM2]=...
                FindOccluPairSurf(defaultPara, ModelInfo1, ModelInfo2, Pair, GlobalScale, Matches, f1, f2)

% Function find the closed occlu point and output a region on epipolar line 
% to search for dense match
tic
% loading Image
Img1 = ModelInfo1.ExifInfo.IDName;
I1 = imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']);
Img2 = ModelInfo2.ExifInfo.IDName;
I2 = imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']);
[Iy1 Ix1] = size(I1);
[DepthY1 DepthX1] = size(ModelInfo1.Model.Depth.FitDepth);
[Iy2 Ix2] = size(I2);
[DepthY2 DepthX2] = size(ModelInfo2.Model.Depth.FitDepth);

% Proper Scaling the Depth 3-d Position and PlaneParameters (all in global scale)
ModelInfo1.Model.PlaneParaInfo.PlanePara = ModelInfo1.Model.PlaneParaInfo.PlanePara/GlobalScale(1);
ModelInfo1.Model.PlaneParaInfo.FitDepth = ModelInfo1.Model.PlaneParaInfo.FitDepth*GlobalScale(1);
ModelInfo1.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo1.Model.PlaneParaInfo.Position3DFited*GlobalScale(1),[ 3 1 2]);
ModelInfo2.Model.PlaneParaInfo.PlanePara = ModelInfo2.Model.PlaneParaInfo.PlanePara/GlobalScale(2);
ModelInfo2.Model.PlaneParaInfo.FitDepth = ModelInfo2.Model.PlaneParaInfo.FitDepth*GlobalScale(1);
ModelInfo2.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo2.Model.PlaneParaInfo.Position3DFited*GlobalScale(2),[ 3 1 2]);
Pair.T = Pair.T/Pair.DepthScale(1)*GlobalScale(1);

% Define possible mathcing point pool
% Specify Target Point in image Pixel position
PointPix1 = f1(:, setdiff(1:length(f1), Matches(1,:))); % unMatched Surf Features
PointPix2 = f2(:, setdiff(1:length(f2), Matches(2,:)));
[Picked_IND1] = ProjPosi2Mask( [Iy1 Ix1], [DepthY1 DepthX1], PointPix1);
[Picked_IND2] = ProjPosi2Mask( [Iy2 Ix2], [DepthY2 DepthX2], PointPix2);

% Region by picking Img1
T1_hat = [[0 -Pair.T(3) Pair.T(2)];...
            [Pair.T(3) 0 -Pair.T(1)];...
            [-Pair.T(2) Pair.T(1) 0]];
F1 = inv(defaultPara.InrinsicK2)'*T1_hat*Pair.R*inv(defaultPara.InrinsicK1);
[Region1 PointPix1 POriReprojM1 PoccluM1]=FindOccluRegion(defaultPara, [Iy2 Ix2], I1, F1, PointPix1, PointPix2, Picked_IND1, Picked_IND2, ModelInfo1, ModelInfo2, Pair);

% Region by picking Img2
Pair2.T = -Pair.R'*Pair.T;
Pair2.R = Pair.R';
Pair2.Xim = Pair.Xim([3 4 1 2],:);
Pair2.DepthScale = Pair.DepthScale([2 1]);
[Region2 PointPix2 POriReprojM2 PoccluM2]=FindOccluRegion(defaultPara, [Iy1 Ix1], I2, F1', PointPix2, PointPix1, Picked_IND2, Picked_IND1, ModelInfo2, ModelInfo1, Pair2);
toc
% Display
if defaultPara.Flag.FindOcclu
        figure;
	MaxD1 = ones(1,size(PoccluM1,2));
	MinD1 = ones(1,size(PoccluM1,2));
	MaxD2 = ones(1,size(PoccluM2,2));
	MinD2 = ones(1,size(PoccluM2,2));
	
	dispMatchSearchRegin(I1, I2,[ PointPix1; ones(1,size(PointPix1,2))], [ PointPix2; ones(1,size(PointPix2,2))], Region1, Region2, F1, ...
        POriReprojM1, MaxD1, PoccluM1, MinD1, ...
        PoccluM2, MaxD2, PoccluM2, MinD2, ...
        0, 'Stacking', 'v', 'Interactive', 0);
end
