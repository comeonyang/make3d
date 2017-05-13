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
function [defaultPara ImgInfo] = PairReInferenceSepRender(defaultPara, Pair, ImgInfo, FlagFirstPair)

% This function is the meta function of 2 decompose inference
% Input:
%	defaultPara - camera intrinsic parameters
%	Pair - Pair image info- matches camera extrinsic
%	ImgInfo - Sup, Ray, MultiScaleSup, SupNeighbor, Depth, Constrain (previous)
% Return:
%	ImgInfo - add new Constrain
% initialize parameter
NegI = diag([1 1 -1]);
NumMatches = length(Pair.Xim);
Img1 = ImgInfo(1).ExifInfo.IDName;
Img2 = ImgInfo(2).ExifInfo.IDName;
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift

% ==========================================================
% ConStrain generation
% 1) Triangulated Depth
x_calib = [ inv(defaultPara.InrinsicK1)*[ Pair.Xim(1:2,:); ones(1, NumMatches)];...
	      inv(defaultPara.InrinsicK2)*[ Pair.Xim(3:4,:); ones(1, NumMatches)]];
TriDepth1 = Pair.lamda(1,:);
TriDepth2 = Pair.lamda(2,:);

% ================= important check scaling stanford ===============
ScaleImg = size(I1);
ScaleDepth = size(ImgInfo(1).Model.Depth.FitDepth);
[IND1] = ProjPosi2Mask( ScaleImg, ScaleDepth, Pair.Xim([1 2],:));
%StichedDeph1 = ImgInfo(1).Model.Depth.RawDepth(IND1);
StichedDeph1 = ImgInfo(1).Model.Depth.FitDepth(IND1);%Min used FitDepth July 1st
if mean(abs(StichedDeph1 - TriDepth1)) < mean( abs( StichedDeph1*Pair.DepthScale(1) - TriDepth1))
    ImgInfo(1).Model.Depth.FitDepth = ImgInfo(1).Model.Depth.FitDepth/Pair.DepthScale(1);
    ImgInfo(1).Model.Depth.RawDepth = ImgInfo(1).Model.Depth.RawDepth/Pair.DepthScale(1);
    disp('Wrong Scaling');
end    

ScaleImg = size(I2);
ScaleDepth = size(ImgInfo(2).Model.Depth.FitDepth);
[IND2] = ProjPosi2Mask( ScaleImg, ScaleDepth, Pair.Xim([3 4],:));
%StichedDeph2 = ImgInfo(2).Model.Depth.RawDepth(IND2);
StichedDeph2 = ImgInfo(2).Model.Depth.FitDepth(IND2);%Min used FitDepth July 1st
if mean(abs(StichedDeph2 - TriDepth2)) < mean( abs( StichedDeph2*Pair.DepthScale(2) - TriDepth2))
    ImgInfo(2).Model.Depth.FitDepth = ImgInfo(2).Model.Depth.FitDepth/Pair.DepthScale(2);
    ImgInfo(2).Model.Depth.RawDepth = ImgInfo(2).Model.Depth.RawDepth/Pair.DepthScale(2);
    disp('Wrong Scaling');    
end 
% ==================================================================
RayNormfactor1 = sqrt( sum(x_calib(1:3,:).^2,1));
RayNormfactor2 = sqrt( sum(x_calib(4:6,:).^2,1));
if FlagFirstPair
	% 2) Ground orientaion and level
	[GroundLevel] = CalGroundLevel(defaultPara, ImgInfo, Pair);
else
    GroundLevel = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'GroundLevel');
end    

[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), Pair.Xim(1:2,:), ImgInfo(1).Model.Depth.FitDepth); % Need IND1
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), Pair.Xim(3:4,:), ImgInfo(1).Model.Depth.FitDepth); % Need IND2
clear I1 I2;

% ==== General Parameters ===============
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
Default.MinTriEffectPercent = 5; % 5 % higher the CoPlaner term
Default.FarestTriDist = 20;%0.1; % unit in pixel; 1 still too samll
%Default.TriCountSupThre = 15;
Default.TriCountSupThre = size(ImgInfo(1).Model.Constrain.RayMatche,1)*0.1;
% ========================================== First Model
    % highly changable
    Default.OutPutFolder = [ defaultPara.OutPutFolder Img1 '/'];
    Default.ScratchFolder = Default.OutPutFolder;
    Default.Wrlname{1} = [defaultPara.Wrlname '_' Img1];
    Default.Flag.AfterInferenceStorage = 0;

    AappendOpt = 0; % not append since generate seperate .wrl files
    ASupMatched = ImgInfo(1).Model.Sup(IND1)';
    mask = ASupMatched == 0;
    ASupMatched(mask)=[];
    ARayMatched = (x_calib(1:3,:)./(repmat( RayNormfactor1,3,1)))';
    ARayMatched(mask,:) = [];
    ADepth_modified = TriDepth1.*RayNormfactor1/Pair.DepthScale(1); % /Pair.DepthScale(1) bing to Ori_scale
    ADepth_modified(:,mask) = [];
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

    Default.filename{1} = [ strrep(ImgInfo(1).ExifInfo.name, '.jpg','') '_']; % Min Modified Aug 18th
%    Default.filename{1} = ['../' strrep(ImgInfo(1).ExifInfo.name, '.jpg','')];
% =====================Meta Model Building =================
if FlagFirstPair
	Rotation = defaultPara.R;
	Translation = defaultPara.T;
	UniScale = defaultPara.Scale;
else
	Rotation = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'R');
	Translation = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'T');
	UniScale = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'Scale');
end

% first Img1
Path = [ defaultPara.OutPutFolder defaultPara.Wrlname '.wrl'];
InLinePath = ['./'  Img1 '/' Default.Wrlname{1} '.wrl'];
R = Rotation;
  % Possibly to be wrong ==========
% 	Angle = recoverAlphasFromU(reshape(R,1,[]));	
% 	Q = GetQauternionFrom2Rotation(zeros(3,1), Angle, false);%[1 0 0 0]';
Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
if any(isnan(Q))
    Q = zeros(4,1);
end    
  % ===============================	
T = Translation;
WRLT = T;
WRLT(3) = -WRLT(3);% for Wrl (-) z component;
Scale = UniScale;
BuildVrmlMetaModel(FlagFirstPair, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
% ============================================================
    if true
     [ ImgInfo(1).Model.PlaneParaInfo] = PlaneParaMRFTriangulateOneShot( Default, Rotation, Translation, AappendOpt, ...
                           [ ARayMatched; Aconstrain.RayMatched],...
                           [ ADepth_modified Aconstrain.Depth_modified]',...
                           [ ASupMatched; Aconstrain.SupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, [], AdepthMap, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, Scale, Rotation, GroundLevel);% eye(3), GroundLevel might be wrong //Min check
    end
    
    % Important Storage the Triangulated info in local scale
    ImgInfo(1).Model.Constrain.RayMatche = [ ARayMatched; Aconstrain.RayMatched];
    ImgInfo(1).Model.Constrain.Depth_modified = [ ADepth_modified Aconstrain.Depth_modified];
    ImgInfo(1).Model.Constrain.SupMatched = [ ASupMatched; Aconstrain.SupMatched];
    model = ImgInfo(1).Model;
    ImgName = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
    save( [defaultPara.ScratchFolder ImgName '/' defaultPara.Wrlname '_' ImgName '_NonMono.mat'], 'model'); % add prefix defaultPara.Wrlname to distinguish model
% =================================== Second Model
    Default.OutPutFolder = [ defaultPara.OutPutFolder Img2 '/'];
    Default.ScratchFolder = Default.OutPutFolder;
    Default.Wrlname{1} = [defaultPara.Wrlname '_' Img2];
    Default.Flag.AfterInferenceStorage = 0;

    Default.RenderFlag =  1;%defaultPara.LastImgFlag;//Min changed to render anyway
    AappendOpt = 0; % not append since generate seperate .wrl files
    ASupMatched = ImgInfo(2).Model.Sup(IND2)';
    mask = ASupMatched == 0;
    ASupMatched(mask)=[];
    ARayMatched = (x_calib(4:6,:)./(repmat( RayNormfactor2,3,1)))';
    ARayMatched(mask,:) = [];
    ADepth_modified = TriDepth2.*RayNormfactor2/Pair.DepthScale(2);% /Pair.DepthScale(2) Bing to Ori_Scale
    ADepth_modified(:,mask) = [];
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

    Default.filename{1} = [ strrep(ImgInfo(2).ExifInfo.name, '.jpg','') '_']; % Min Modified Aug 18th
%    Default.filename{1} = ['../' strrep(ImgInfo(2).ExifInfo.name, '.jpg','')];
% ========Add model or modify model in Global Meta Model file=================
% Second Img2
Path = [ defaultPara.OutPutFolder defaultPara.Wrlname '.wrl'];
InLinePath = ['./'  Img2 '/' Default.Wrlname{1} '.wrl'];
R = Rotation*Pair.R';
Q = Rotation2Q(NegI*R*NegI);
if any(isnan(Q))
    Q = zeros(4,1);
end 
  % ===============================	
T = Translation + Rotation*(-Pair.R'*(Pair.T/Pair.DepthScale(1)*UniScale ));
WRLT = T;
WRLT(3) = -WRLT(3);% for Wrl (-) z component;
Scale = Pair.DepthScale(2)/Pair.DepthScale(1)*UniScale;
BuildVrmlMetaModel(0, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
% ============================================================================true
if true;%Default.RenderFlag % Since using inline to link VRML render it anyway
     [ ImgInfo(2).Model.PlaneParaInfo] = PlaneParaMRFTriangulateOneShot( Default, R, T, AappendOpt, ...
                           [ ARayMatched; Aconstrain.RayMatched],...
                           [ ADepth_modified Aconstrain.Depth_modified]',...
                           [ ASupMatched; Aconstrain.SupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, [], AdepthMap, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, Scale, R, GroundLevel); % eye(3), GroundLevel might be wrong //Min check
end

    % Important Storage the Triangulated info in local scale
    ImgInfo(2).Model.Constrain.RayMatche = [ ARayMatched; Aconstrain.RayMatched];
    ImgInfo(2).Model.Constrain.Depth_modified = [ ADepth_modified Aconstrain.Depth_modified];
    ImgInfo(2).Model.Constrain.SupMatched = [ ASupMatched; Aconstrain.SupMatched];
    model = ImgInfo(2).Model;
    ImgName = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
    save( [defaultPara.ScratchFolder ImgName '/' defaultPara.Wrlname '_' ImgName '_NonMono.mat'], 'model');% add prefix defaultPara.Wrlname to distinguish model


% Save GroundLevel, R, T, and, Scale info in /data/ModelStatus for adding
% new  image into the current model 
SaveModelSatus( defaultPara, defaultPara.Wrlname, {Img1,Img2}, GroundLevel, R, T, Scale, FlagFirstPair)

return;

