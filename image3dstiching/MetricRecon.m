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
function [defaultPara ImgInfo fail] = MetricRecon(defaultPara, ImgInfo)

% This function do Metric Reconstruction for a pair of images
% Step outline:
% 	1) Mono calulation
%	2) Robust Pose and Match generation
%	3) Dense match 
%	4) Re-inference
%	5) Coherence Texture Process 
% 	6) Volumatric refining
%	7) Rendering
%
% Input:
%	1) defaultPara - camera intrinsic, ....etc.
%	2) ImgInfo - special structure with Exif, GPS, IMU, Model info
%
% Return:
%	1) TriPoint - Triangluated point

Debug = 0;
NumPcik = 10;
fail= 0;

% Parameter Setting
Img1 = ImgInfo(1).ExifInfo.IDName;
Img2 = ImgInfo(2).ExifInfo.IDName;
disp([ 'Processing ' Img1 ' & ' Img2]);
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift		
ImgScale1 = size(I1);
ImgScale2 = size(I2);
[Img1Index] = ImgInfoIndexFromName(ImgInfo, Img1);
[Img2Index] = ImgInfoIndexFromName(ImgInfo, Img2);
NegI = diag([1 1 -1]);

% 1. Mono calulation or load the pre-calculated data ------------------------
[ ImgInfo] = SingleModelInfo(defaultPara, ImgInfo);

% 2. Robust Pose and Match generation ---------------------------------------------
if defaultPara.Flag.PoseMatches
	[Pair, ImgInfo, Matches, fail] = PoseMatchEst(defaultPara, ImgInfo);
	Pair.matches = Matches;
	if fail > 0		
		disp('End of Metric Reconstruction');
		return;
	end
	if defaultPara.Flag.PoseMatchStor
		save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_PoseMatch.mat'],'Pair','Matches','fail');
	end
end
% To Check if Placing MonoModel using R T make sense
% Indep_Mono_Model_CheckRT(defaultPara, ImgInfo, R, T);

% 3. Re-inference -----------------------------------------------------------------
[status1, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '/' defaultPara.Wrlname '_' Img1 '_NonMono.mat']);
[status2, result] = system(['ls ' defaultPara.Fdir '/data/' Img2 '/' defaultPara.Wrlname '_' Img2 '_NonMono.mat']);
if defaultPara.Flag.ReInference || status1 || status2
	if ~defaultPara.Flag.PoseMatches
		[Pair, Matches, fail] = LoadPoseMatch(defaultPara.Fdir, Img1, Img2);
		Pair.matches = Matches;
	end
     [defaultPara ImgInfo] = PairReInferenceSepRender(defaultPara, Pair, ImgInfo, defaultPara.Flag.FlagFirstPair);
end

% 4. Pair-Wise Occlusion detection and Corrolation Matching ----------------------------------------------------
if defaultPara.Flag.Refinement
	if ~defaultPara.Flag.ReInference && ~defaultPara.Flag.PoseMatches
		[ImgInfo1 ImgInfo2 Img1Index Img2Index Pair GlobalScale] = ... % only deal with pairs have been matched
			LoadDataForFindOcclu(defaultPara, defaultPara.Wrlname, Img1, Img2, ...
			defaultPara.Flag.FlagImgInfoLoadPreStorage); % load all imformation needed for occlusion detection
		ImgInfo = [ImgInfo1 ImgInfo2];
	else
		GlobalScale(1) = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'Scale');
		GlobalScale(2) = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img2, 'Scale');
	end

	% First to Occlusion detection
	% time consuming about 5mins
	disp('Start Surf Feature Point Occlusion detection');
	[status, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Match.mat']);
	if ~defaultPara.Flag.FlagPreloadOccluDetectMatches || status
		s = warning('off');
		[PointPix1 PointDepth1 FaceSetPickedIND1 POriReprojM1 FieldOccluPix1 OccluDist1 OccluedFaceSetIND1 OccluedFaceSetIDRemained1...
		PointPix2 PointDepth2 FaceSetPickedIND2 POriReprojM2 FieldOccluPix2 OccluDist2 OccluedFaceSetIND2 OccluedFaceSetIDRemained2] = ...
		FindOccluPair(defaultPara, ImgInfo(1), ImgInfo(2), Pair, GlobalScale, 1); % detect occlsion (use Ray as detector not the surf Features)
		% data Struction define:
		% PointPix PointDepth FaceSetPickedIND POriReprojM FieldOccluPix OccluDist OccluedFaceSetIND
		% -- allof the size that single ray pass through both FaceSet
		% FaceSetPickedIND1 used in DepthMap size or surfFeature size
		% OccluedFaceSetIND1 used in DepthMap size
		save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Match.mat'],...
		'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1',...
			'OccluedFaceSetIND1','OccluedFaceSetIDRemained1', ...
		'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2',...
			'OccluedFaceSetIND2','OccluedFaceSetIDRemained2');
		s = warning('on');
	else
		load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Match.mat']);
	end       
	
	% Define occlusion if OccluDist > CentainThreshold
	if defaultPara.Flag.FlagRecipicalOccluDetection
		Mask1 = PointDepth1./OccluDist1 > defaultPara.OcclusionRatioThre | OccluDist1./PointDepth1 > defaultPara.OcclusionRatioThre;
		Mask2 = PointDepth2./OccluDist2 > defaultPara.OcclusionRatioThre | OccluDist2./PointDepth2 > defaultPara.OcclusionRatioThre;
	else
		Mask1 = PointDepth1./OccluDist1 > defaultPara.OcclusionRatioThre;%defaultPara.OccluDistThre;
		Mask2 = PointDepth2./OccluDist2 > defaultPara.OcclusionRatioThre;%defaultPara.OccluDistThre;
	end

	% Second do Surf Features Matches given Occlusion infomation 
	disp('Start Finding Surf Feature Point Occlusion Match');
	[status, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '_' Img2 '_SurfOccluMatches.mat']);
    if defaultPara.Flag.LoadStoragedSurfOccluMatches && ~status
        load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_SurfOccluMatches.mat']);
        SurfMatches = matches;
    else
        [SurfMatches fail] = MatchPointsGivenOcclusion(defaultPara, ImgInfo(1), ImgInfo(2), ImgScale1, ImgScale2, ...
            Img1, Img2, Img1Index, Img2Index, ...
            Pair,  GlobalScale, defaultPara.Wrlname, 'SurfOccluMatches', ...
            POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), FaceSetPickedIND1(:,Mask1), ...
            POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), FaceSetPickedIND2(:,Mask2), true );	
    end    
    if fail > 0		
		disp('End of Metric Reconstruction');
		return;
    end
    
	% 6. Pair-Wise Robust Pose and Match generation (Second time) given the new Occlusion Matches
	% ===== Debug Only
	if defaultPara.Flag.FlagStorageBeforePairNew
		disp('Storaging Data Before Runing PairNew');
		save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_BeforePairNew.mat']);
	end
	% ===============
	PriorPose.R = [Pair.R; Pair.R'];
	PriorPose.T = [Pair.T; -Pair.R'*Pair.T];
    [status, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '_' Img2 '_PairNew.mat']);
    if ~defaultPara.Flag.loadStoragedPairNew || status
        [PairNew fail]=TestWholePostMatch(defaultPara, ImgInfo(1), ImgInfo(2), PriorPose);
        if defaultPara.Flag.FlagStorageAfterPairNew           
            save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_PairNew.mat'],'PairNew','fail');
        end
    else
        load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_PairNew.mat']);
    end    
    if fail > 0		
		disp('End of Metric Reconstruction');
		return;
    end
    
	% Given PairNew info do Correlation Matches
	% Need to detect Occlusion for Correlation matches
	% time consuming about 5mins
	disp('Start Ray Point Occlusion detection');
	[status, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Ray.mat']);
	if ~defaultPara.Flag.FlagPreloadOccluDetectRay || status
		s = warning('off');
		[PointPix1 PointDepth1 FaceSetPickedIND1 POriReprojM1 FieldOccluPix1 OccluDist1 OccluedFaceSetIND1 OccluedFaceSetIDRemained1...
		PointPix2 PointDepth2 FaceSetPickedIND2 POriReprojM2 FieldOccluPix2 OccluDist2 OccluedFaceSetIND2 OccluedFaceSetIDRemained2] = ...
		FindOccluPair(defaultPara, ImgInfo(1), ImgInfo(2), Pair, GlobalScale, false); % detect occlsion (use Ray as detector not the surf Features)
		% data Struction define:
		% PointPix PointDepth FaceSetPickedIND POriReprojM FieldOccluPix OccluDist OccluedFaceSetIND
		% -- allof the size that single ray pass through both FaceSet
		% FaceSetPickedIND1 used in DepthMap size or surfFeature size
		% OccluedFaceSetIND1 used in DepthMap size
		save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Ray.mat'],...
		'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1',...
			'OccluedFaceSetIND1','OccluedFaceSetIDRemained1', ...
		'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2',...
			'OccluedFaceSetIND2','OccluedFaceSetIDRemained2');
		s = warning('on');
	else
		load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' defaultPara.PostFixStrAfter '_OccluDetect_Ray.mat']);
	end       
		
	% Define occlusion if OccluDist > CentainThreshold
	if defaultPara.Flag.FlagRecipicalOccluDetection
		Mask1 = PointDepth1./OccluDist1 > defaultPara.OcclusionRatioThre | OccluDist1./PointDepth1 > defaultPara.OcclusionRatioThre;
		Mask2 = PointDepth2./OccluDist2 > defaultPara.OcclusionRatioThre | OccluDist2./PointDepth2 > defaultPara.OcclusionRatioThre;
	else
		Mask1 = PointDepth1./OccluDist1 > defaultPara.OcclusionRatioThre;%defaultPara.OccluDistThre;
		Mask2 = PointDepth2./OccluDist2 > defaultPara.OcclusionRatioThre;%defaultPara.OccluDistThre;
	end

	disp('Start Finding Ray Point Occlusion Match');
	[status, result] = system(['ls ' defaultPara.Fdir '/data/' Img1 '_' Img2 '_CorrMatches.mat']);
	if defaultPara.Flag.LoadCorrMatches && ~status
		load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_CorrMatches.mat']);
		CorrMatches = Matches;
    else
        % debuging purpose ==========
        if Debug
            Mask1 = find(Mask1);
            Mask2 = find(Mask2);            
            RandPcik1 = randperm(length(Mask1));
            NumPcik1 = min(NumPcik,length(Mask1));
            Mask1 = Mask1(RandPcik1(1:NumPcik1));
            RandPcik2 = randperm(length(Mask2));
            NumPcik2 = min(NumPcik,length(Mask2));
            Mask2 = Mask2(RandPcik2(1:NumPcik));
        end    
        % ===========================
		CorrTime = tic;
		fprintf(' Running CorrMatches ........');
		[ CorrMatches]=CorrMatchPointsGivenOcclusion( ...
			defaultPara, PairNew, GlobalScale,...
			POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), PointPix1(:,Mask1),...
			POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), PointPix2(:,Mask2),...
			ImgInfo(1), ImgInfo(2), ImgScale1, ImgScale2, ...
			true);
		disp(['		' num2str( toc( CorrTime)) ' seconds.']);
	end

	% Remove redundant matches in SurfMatches & CorrMatches
	Coeff = [ones(1, size(SurfMatches, 2)) zeros(1, size(CorrMatches, 2))];
	AllMatches = [ SurfMatches CorrMatches];
	[Inliers] = CleanMatch(AllMatches, Coeff);
	[ReverseInliers] = CleanMatch(AllMatches([3 4 1 2], Inliers), Coeff(Inliers));
	Inliers = Inliers(ReverseInliers);
	AllMatches = AllMatches(:,Inliers);	

	% triangulate depth and put into constrain for inference
	if ~isempty(AllMatches)%
		disp('Start Triangluation and Storage constrain for inference');
		PostFixStrAfter = 'NewPair';
		tempf1 = AllMatches(1:2,:);
		tempf2 = AllMatches(3:4,:);
		x_calib = [ inv(defaultPara.InrinsicK1)*[ tempf1; ones(1,size(tempf1,2))];...
		inv(defaultPara.InrinsicK2)*[ tempf2; ones(1,size(tempf2,2))]];
		[ lamda1 lamda2 Error] = triangulation( defaultPara, Pair.R, Pair.T, x_calib);
		% notice lamda re-scale to local model scale
		lamda1 = lamda1./GlobalScale(1);
		lamda2 = lamda2./GlobalScale(2);
		% Storage the match result for later ReInference
		AddMatch2Model(defaultPara, defaultPara.Wrlname, lamda1, AllMatches(1:2,:), ImgInfo(1), ImgScale1, Img1Index, Img2Index, PostFixStrAfter, Error);
		AddMatch2Model(defaultPara, defaultPara.Wrlname, lamda2, AllMatches(3:4,:), ImgInfo(2), ImgScale2, Img2Index, Img1Index, PostFixStrAfter, Error);
		% Storage the New Matches                   
		if defaultPara.Flag.FlagRefinementDisp
			disp('Storaging Occlusion Correlation Matches');
			save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_AllMatches.mat'],'AllMatches','Error');
		end
	end


	% 7. Re-inference ------ with new Pose and Matchese
	if defaultPara.Flag.RefinementInference
		% For Img1 ======================================================
		disp('PaieNew Inference for Img1');
		% Data preparing and ReInfernece
		Default.Wrlname{1} = [defaultPara.Wrlname '_' Img1 '_' PostFixStrAfter];   
		load([defaultPara.Fdir '/data/' Img1 '/' defaultPara.Wrlname '_' Img1 '_' PostFixStrAfter '.mat']);
		R = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'R');
		T = LoadModelStatus( defaultPara.Fdir, defaultPara.Wrlname, Img1, 'T');
		Scale = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'Scale');
		GroundLevel = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'GroundLevel');
		SingleImgReInference(defaultPara, model, Img1, GroundLevel, defaultPara.Wrlname, ...
			R, T, Scale ,PostFixStrAfter);

		%   Build Meta Wrl file =========
		Path = [ defaultPara.OutPutFolder defaultPara.Wrlname '_' PostFixStrAfter '.wrl'];
		InLinePath = ['./'  Img1 '/' Default.Wrlname{1} '.wrl'];
		Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
		if any(isnan(Q))
			Q = zeros(4,1);
		end
		WRLT = T;
		WRLT(3) = -WRLT(3);% for Wrl (-) z component;
		Scale = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img1, 'Scale');
 		[status, result] = system(['ls ' Path]);
		if status
			BuildVrmlMetaModel(1, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
		else
			BuildVrmlMetaModel(0, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
		end
		% ===============================

		% For Img2 ===================================
		disp('PaieNew Inference for Img2');
		% Data preparing and ReInfernece
		Default.Wrlname{1} = [defaultPara.Wrlname '_' Img2 '_' PostFixStrAfter];   
		load([defaultPara.Fdir '/data/' Img2 '/' defaultPara.Wrlname '_' Img2 '_' PostFixStrAfter '.mat']);
		T = T + R*(-PairNew.R'*(PairNew.T/PairNew.DepthScale(1)*Scale ));
		R = R*PairNew.R';
		Scale = PairNew.DepthScale(2)/PairNew.DepthScale(1)*Scale;
		SingleImgReInference(defaultPara, model, Img2, GroundLevel, defaultPara.Wrlname, ...
			R, T, Scale ,PostFixStrAfter);

		%   Build Meta Wrl file =========
		Path = [ defaultPara.OutPutFolder defaultPara.Wrlname '_' PostFixStrAfter '.wrl'];
		InLinePath = ['./'  Img2 '/' Default.Wrlname{1} '.wrl'];
		Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
		if any(isnan(Q))
			Q = zeros(4,1);
		end
		WRLT = T;
		WRLT(3) = -WRLT(3);% for Wrl (-) z component;
		Scale = LoadModelStatus(defaultPara.Fdir, defaultPara.Wrlname, Img2, 'Scale');
 		[status, result] = system(['ls ' Path]);
		if status
			BuildVrmlMetaModel(1, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
		else
			BuildVrmlMetaModel(0, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
		end
		% Storge the New GlobalPose info of Img2
		SaveModelSatus( defaultPara, defaultPara.Wrlname, {Img1,Img2}, GroundLevel, R, T, Scale, defaultPara.Flag.FlagFirstPair);
		% ==============================
	end
end

disp('End of Metric Reconstruction');
return;

