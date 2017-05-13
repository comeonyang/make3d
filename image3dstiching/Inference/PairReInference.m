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
function [defaultPara ImgInfo GlabalInfo] = PairReInference(defaultPara, Pair, ImgInfo, FlagFirstPair)

% This function is the meta function of 2 decompose inference
% Input:
%	defaultPara - camera intrinsic parameters
%	Pair - Pair image info- matches camera extrinsic
%	ImgInfo - Sup, Ray, MultiScaleSup, SupNeighbor, Depth, Constrain (previous)
% Return:
%	ImgInfo - add new Constrain
% initialize parameter
NumMatches = length(Pair.Xim);
Img1 = ImgInfo(1).ExifInfo.IDName;
Img2 = ImgInfo(2).ExifInfo.IDName;
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift

% ===================== Rescale Depth and T prorperly ======
if FlagFirstPair
	UniScale = defaultPara.Scale;
else
	UniScale = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'Scale');
end
ImgInfo(1).Model.Depth.FitDepth = ImgInfo(1).Model.Depth.FitDepth*UniScale;
ImgInfo(1).Model.Depth.RawDepth = ImgInfo(1).Model.Depth.RawDepth*UniScale;
ImgInfo(2).Model.Depth.FitDepth = ImgInfo(2).Model.Depth.FitDepth*Pair.DepthScale(2)/Pair.DepthScale(1)*UniScale;
ImgInfo(2).Model.Depth.RawDepth = ImgInfo(2).Model.Depth.RawDepth*Pair.DepthScale(2)/Pair.DepthScale(1)*UniScale;
Pair.T = Pair.T/Pair.DepthScale(1)*UniScale;
Scale = Pair.DepthScale(2)/Pair.DepthScale(1)*UniScale;

% ==========================================================
% ConStrain generation
% 1) Triangulated Depth
x_calib = [ inv(defaultPara.InrinsicK1)*[ Pair.Xim(1:2,:); ones(1, NumMatches)];...
	      inv(defaultPara.InrinsicK2)*[ Pair.Xim(3:4,:); ones(1, NumMatches)]];
[ TriDepth1 TriDepth2] = triangulation( defaultPara, Pair.R, Pair.T, x_calib);
%TriDepth1 = Pair.lamda(1);
%TriDepth2 = Pair.lamda(2);
RayNormfactor1 = sqrt( sum(x_calib(1:3,:).^2,1));
RayNormfactor2 = sqrt( sum(x_calib(4:6,:).^2,1));
if FlagFirstPair
	% 2) Ground orientaion and level
	[GroundLevel] = CalGroundLevel(defaultPara, ImgInfo, Pair);
else
    GroundLevel = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'GroundLevel');
end    

% ImgInfo(1) ReInference
%[ImgInfo(1)] = ReInference(defaultPara, Pair.R', -R'*Pair.T, Pair.Xim(1:4, :), ImgInfo([1 2]), [ TriDepth1 TriDepth2], GroundLevel, eye(3), 1);
% ImgInfo(2) ReInference
%[ImgInfo(2)] = ReInference(defaultPara, Pair.R, Pair.T, Pair.Xim([3 4 1 2], :), ImgInfo([2 1]), [ TriDepth2 TriDepth1], GroundLevel, eye(3), FlagRenderingAll);

% ================================old method
[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), Pair.Xim(1:2,:), ImgInfo(1).Model.Depth.FitDepth);
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), Pair.Xim(3:4,:), ImgInfo(1).Model.Depth.FitDepth);
clear I1 I2;
    Default.OutPutFolder = defaultPara.OutPutFolder;
    Default.ScratchFolder = defaultPara.ScratchFolder;
    Default.Wrlname{1} = defaultPara.Wrlname;
    Default.Flag.AfterInferenceStorage = 0;
Default.RenderFlag =  defaultPara.RenderFlag;
[Default.VertYNuDepth Default.HoriXNuDepth] = size(ImgInfo(1).Model.Depth.FitDepth);
Default.fy = 2400.2091651084;
Default.fx = 2407.3312729885838;
Default.Ox = 1110.7122391785729;%2272/2; %
Default.Oy = 833.72104535435108;%1704/2; %
Default.a_default = 2272/Default.fx; %0.70783777; %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
Default.b_default = 1704/Default.fy; %0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
Default.Ox_default = 1-Default.Ox/2272;%0.489272914; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;%0.488886982; % camera origin offset from the image center in vertical direction
Default.MinTriEffectPercent = 5;
Default.FarestTriDist = 0.1;
Default.TriCountSupThre = 10;

    AappendOpt = ~FlagFirstPair;
    ASupMatched = ImgInfo(1).Model.Sup(IND1)';
    mask = ASupMatched == 0;
    ASupMatched(mask)=[];
    ARayMatched = (x_calib(1:3,:)./(repmat( RayNormfactor1,3,1)))';
    ARayMatched(mask,:) = [];
ADepth_modified = TriDepth1.*RayNormfactor1;	
    ADepth_modified(:,mask) = [];
if FlagFirstPair
	ARotation = defaultPara.R;
	ATranslation = defaultPara.T;
else
	ARotation = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'R');
	ATranslation = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'T');
end
Aconstrain.RayMatched = ImgInfo(1).Model.Constrain.RayMatche;
Aconstrain.Depth_modified = ImgInfo(1).Model.Constrain.Depth_modified;
Aconstrain.SupMatched = ImgInfo(1).Model.Constrain.SupMatched;
ASup = ImgInfo(1).Model.Sup;
ASupOri = ASup;
AdepthMap = ImgInfo(1).Model.Depth.RawDepth; % or RawDepth
ARayOri = ImgInfo(1).Model.Ray;
ARayAll = ARayOri;
ASupNeighborTable = ImgInfo(1).Model.SupNeighborTable;
AmaskSky = ImgInfo(1).Model.maskSky;
AmaskG = ImgInfo(1).Model.maskG;
AMultiScaleSupTable = ImgInfo(1).Model.MultiScaleSupTable;

    Default.filename{1} = strrep(ImgInfo(1).ExifInfo.name, '.jpg','');
if true
     [ ImgInfo(1).Model.PlaneParaInfo] = PlaneParaMRFTriangulateOneShot( Default, ARotation, ATranslation, AappendOpt, ...
                           [ ARayMatched; Aconstrain.RayMatched],...
                           [ ADepth_modified Aconstrain.Depth_modified]',...
                           [ ASupMatched; Aconstrain.SupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, [], AdepthMap, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, 0, ARotation, GroundLevel);
end
ImgInfo(1).Model.Constrain.RayMatche = [ ARayMatched; Aconstrain.RayMatched];
ImgInfo(1).Model.Constrain.Depth_modified = [ ADepth_modified Aconstrain.Depth_modified];
ImgInfo(1).Model.Constrain.SupMatched = [ ASupMatched; Aconstrain.SupMatched];
model = ImgInfo(1).Model;
% model.PlaneParaInfo = ImgInfo(1).PlaneParaInfo;
ImgName = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
save( [defaultPara.ScratchFolder ImgName '/' ImgName '_NonMono.mat'], 'model');
GlabalInfo(2).GroundLevel = GroundLevel;
GlabalInfo(2).ARotation = ARotation;
GlabalInfo(2).ATranslation = ATranslation;
GlabalInfo(2).Scale = UniScale;

% ===================================
Default.RenderFlag =  defaultPara.LastImgFlag;
    AappendOpt = 1;
    ASupMatched = ImgInfo(2).Model.Sup(IND2)';
    mask = ASupMatched == 0;
    ASupMatched(mask)=[];
    ARayMatched = (x_calib(4:6,:)./(repmat( RayNormfactor2,3,1)))';
    ARayMatched(mask,:) = [];
ADepth_modified = TriDepth2.*RayNormfactor2;	
    ADepth_modified(:,mask) = [];
ARotation1 = ARotation; % Bad method to avoid override
ATranslation1 = ATranslation; % Bad method to avoid override
ARotation = ARotation1*Pair.R';
ATranslation = ATranslation1 + ARotation1*(-Pair.R'*Pair.T);
Aconstrain.RayMatched = ImgInfo(2).Model.Constrain.RayMatche;
Aconstrain.Depth_modified = ImgInfo(2).Model.Constrain.Depth_modified;
Aconstrain.SupMatched = ImgInfo(2).Model.Constrain.SupMatched;
ASup = ImgInfo(2).Model.Sup;
ASupOri = ASup;
AdepthMap = ImgInfo(2).Model.Depth.RawDepth; % or RawDepth
ARayOri = ImgInfo(2).Model.Ray;
ARayAll = ARayOri;
ASupNeighborTable = ImgInfo(2).Model.SupNeighborTable;
AmaskSky = ImgInfo(2).Model.maskSky;
AmaskG = ImgInfo(2).Model.maskG;
AMultiScaleSupTable = ImgInfo(2).Model.MultiScaleSupTable;

    Default.filename{1} = strrep(ImgInfo(2).ExifInfo.name, '.jpg','');
if Default.RenderFlag
     [ ImgInfo(2).Model.PlaneParaInfo] = PlaneParaMRFTriangulateOneShot( Default, ARotation, ATranslation, AappendOpt, ...
                           [ ARayMatched; Aconstrain.RayMatched],...
                           [ ADepth_modified Aconstrain.Depth_modified]',...
                           [ ASupMatched; Aconstrain.SupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, [], AdepthMap, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, 0, ARotation, GroundLevel);
end

ImgInfo(2).Model.Constrain.RayMatche = [ ARayMatched; Aconstrain.RayMatched];
ImgInfo(2).Model.Constrain.Depth_modified = [ ADepth_modified Aconstrain.Depth_modified];
ImgInfo(2).Model.Constrain.SupMatched = [ ASupMatched; Aconstrain.SupMatched];
model = ImgInfo(2).Model;
% model.PlaneParaInfo = ImgInfo(2).PlaneParaInfo;
ImgName = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
save( [defaultPara.ScratchFolder ImgName '/' ImgName '_NonMono.mat'], 'model');

% modify defaultPara
%defaultPara.GroundLevel = GroundLevel;
%defaultPara.R = ARotation;%defaultPara.R*Pair.R';
%defaultPara.T = ATranslation;%defaultPara.T + defaultPara.R*(-Pair.R'*Pair.T);
%defaultPara.Scale = Pair.DepthScale(2)/Pair.DepthScale(1)*UniScale;
SaveModelSatus( defaultPara.Fdir, defaultPara.Wrlname, Img2, GroundLevel, ARotation, ATranslation, Scale, FlagFirstPair)
GlabalInfo(2).GroundLevel = GroundLevel;
GlabalInfo(2).ARotation = ARotation;
GlabalInfo(2).ATranslation = ATranslation;
GlabalInfo(2).Scale = Scale;

return;

