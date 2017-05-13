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
function [] = OcclusionRemoveRefineMent(defaultPara, Wrlname, PairList, PostFixStrAfter, Type)

% This function do three thing
% 1) detect the occlusion from the current model
% 2) decide redundent plane
% 3) Re-Render by not rendering the redundent plane information

% Input:
% 1) defaultPara - all parameter share for the whole stitching3d folder
%	2) Wrlname - the name of the model, different by the PairList
%	3) PairList - the linearly order of image that been added to the model

% Fixed constant --- need to be put in defaultPara by Min later
NegI = diag([1 1 -1]);
PostFixStr = 'NonMono';%'NonMonoOccluMatched'
if nargin <4
    PostFixStrAfter = 'OccluRemove';
end    
DisCtsThre = 0.1; % fractional error
OcclusionRatioThre = 1.1;
SurfFlag = 0; % default not using surf Feature point to find occlusion
if strcmp( lower(Type), 'matches')
	SurfFlag = 1;
end

% 1) Run occlusion detection and Matcheing for every pair of images

ImgList = unique(PairList);
NumImg = length(ImgList);

% copy the ImgInfo of each image to another "PostFixStrAfter".mat file
for i = 1:NumImg
	load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStr '.mat']);
	save([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat'],'model');	
end

%  	detect occlusion for every pair of image
for i = 1:NumImg
	for j = (i+1):NumImg

%  	[Img1Index] = ImgInfoIndexFromName(ImgInfo, ImgList{i});
%  	[Img2Index] = ImgInfoIndexFromName(ImgInfo, ImgList{j});
	
	[ImgInfo1 ImgInfo2 Img1Index Img2Index Pair GlobalScale] = ... % What to do with Pairs haven't been matched ???????????????????
		LoadDataForFindOcclu(defaultPara, Wrlname, ImgList{i}, ImgList{j}, ...
                		     0, PostFixStrAfter); % load all imformation needed for occlusion detection

    % Define variables
	H = 2274;
        V = 1704;
        Img1 = ImgInfo1.ExifInfo.IDName;
        Img2 = ImgInfo2.ExifInfo.IDName;
        I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
        I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift		
	ImgScale1 = size(I1);
        ImgScale2 = size(I2);
    
    % time consuming about 5mins
    if ~defaultPara.Flag.FlagPreloadOccluDetect
        [PointPix1 PointDepth1 FaceSetPickedIND1 POriReprojM1 FieldOccluPix1 OccluDist1 OccluedFaceSetIND1 ...
            PointPix2 PointDepth2 FaceSetPickedIND2 POriReprojM2 FieldOccluPix2 OccluDist2 OccluedFaceSetIND2] = ...
            FindOccluPair(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, SurfFlag); % detect occlsion (use Ray as detector not the surf Features)
	% data Struction define:
	% PointPix PointDepth FaceSetPickedIND POriReprojM FieldOccluPix OccluDist OccluedFaceSetIND
	% -- allof the size that single ray pass through both FaceSet
	% FaceSetPickedIND1 used in DepthMap size or surfFeature size
 	% OccluedFaceSetIND1 used in DepthMap size
        if SurfFlag
            save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect_Match.mat'],...
                'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1','OccluedFaceSetIND1',...
                'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2','OccluedFaceSetIND2');
        else    
            save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect_Ray.mat'],...
                'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1','OccluedFaceSetIND1',...
                'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2','OccluedFaceSetIND2');
        end    
    else
        if SurfFlag
            load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect_Match.mat']);
        else    
            load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect_Ray.mat']);
        end    
    end        

	% Define occlusion if OccluDist > CentainThreshold
	Mask1 = PointDepth1./OccluDist1 > OcclusionRatioThre;%defaultPara.OccluDistThre;
	Mask2 = PointDepth2./OccluDist2 > OcclusionRatioThre;%defaultPara.OccluDistThre;
if false
	% Action: find new matches or decide plane to remove given the occlusion infomation
	switch lower(Type)
	
	case 'matches'
		MatchPointsGivenOcclusion(defaultPara, ImgScale1, ImgScale2, Img1, Img2, Img1Index, Img2Index, Pair,  ...
			POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), FaceSetPickedIND1(:,Mask1), ...
			POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), FaceSetPickedIND2(:,Mask2) );	
	case 'Remove'
		RemovePlaneGivenOcclusion(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, ...
			FaceSetPickedIND1, FaceSetPickedIND1, FaceSetPickedIND2, OccluedFaceSetIND2, Mask1, Mask2);
    end
end

    end
	
end

% 2) Render of each image individiually

for i = 1:NumImg
	
%	Data preparing and ReInfernece
    Default.Wrlname{1} = [Wrlname '_' ImgList{i} '_' PostFixStrAfter];   
    Default.filename{1} = ['../' ImgList{i}];
    load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat']);
    Default.OutPutFolder = [ defaultPara.OutPutFolder ImgList{i} '/'];
        
    if false
        DiffHori = conv2(model.PlaneParaInfo.FitDepth, [1 -1], 'valid');
        DiffVert = conv2(model.PlaneParaInfo.FitDepth, [1; -1], 'valid');
        maskH = [ (DiffHori./model.PlaneParaInfo.FitDepth(:,2:end) > DisCtsThre) zeros(size(model.PlaneParaInfo.FitDepth,1),1)];
        maskH(end,:) = 0;
        maskV = [ (DiffVert./model.PlaneParaInfo.FitDepth(2:end,:) > DisCtsThre); zeros(1, size(model.PlaneParaInfo.FitDepth,2))];
        maskV(:,end) = 0;
        model.PlaneParaInfo.SupOri(logical(  maskH) ) = 0;
        model.PlaneParaInfo.SupOri(logical(  maskV) ) = 0;
    end

	WrlFacestHroiReduce(model.PlaneParaInfo.Position3DFited, model.PlaneParaInfo.PositionTex, model.PlaneParaInfo.SupOri, ...
		Default.filename{1}, Default.Wrlname{1}, ...
        	Default.OutPutFolder, 0, 0);
%	SingleImgReInference(defaultPara, model, ImgList{i}, LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'GroundLevel'));
    
%   Build Meta Wrl file

    Path = [ defaultPara.OutPutFolder Wrlname '_' PostFixStrAfter '.wrl'];
    InLinePath = ['./'  ImgList{i} '/' Default.Wrlname{1} '.wrl'];
    R = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'R');
    Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
    if any(isnan(Q))
        Q = zeros(4,1);
    end
    % ===============================
    T = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'T');
    WRLT = T;
    WRLT(3) = -WRLT(3);% for Wrl (-) z component;
    Scale = LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'Scale');
    if i == 1
        BuildVrmlMetaModel(1, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
    else
        BuildVrmlMetaModel(0, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1)); 
    end
end

return;
