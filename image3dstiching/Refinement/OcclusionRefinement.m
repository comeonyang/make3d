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
function [] = OcclusionRefinement(defaultPara, Wrlname, PairList, PostFixStrAfter, Type, PostFixStr)

% This function do three thing
% 1) detect the occlusion from the current model
% 2) decide redundent plane
% 3) Re-Render by not rendering the redundent plane information

% Input:
% 1) defaultPara - all parameter share for the whole stitching3d folder
%	2) Wrlname - the name of the model, different by the PairList
%	3) PairList - the linearly order of image that been added to the model

% Fixed constant --- need to be put in defaultPara by Min later
depthratioMin = 0.01;
depthratioMax = 100;
NegI = diag([1 1 -1]);
if nargin <6
    PostFixStr = 'NonMono';%'NonMonoOccluMatched'
end    
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

if defaultPara.Flag.FlagOccluDetect 

	tic; % start timing the time for the OccluDetection Step

	% copy the ImgInfo of each image to another "PostFixStrAfter".mat file
    	if defaultPara.Flag.FlagResetModeltoNomoModel
        	for i = 1:NumImg
            		load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStr '.mat']);
            		save([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat'],'model');	
        	end
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
	        	Img1 = ImgInfo1.ExifInfo.IDName
		        Img2 = ImgInfo2.ExifInfo.IDName
        		I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
		        I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift		
			ImgScale1 = size(I1);
	        	ImgScale2 = size(I2);
    
	    		% time consuming about 5mins
		    	if ~defaultPara.Flag.FlagPreloadOccluDetect
        			[PointPix1 PointDepth1 FaceSetPickedIND1 POriReprojM1 FieldOccluPix1 OccluDist1 OccluedFaceSetIND1 OccluedFaceSetIDRemained1...
            			PointPix2 PointDepth2 FaceSetPickedIND2 POriReprojM2 FieldOccluPix2 OccluDist2 OccluedFaceSetIND2 OccluedFaceSetIDRemained2] = ...
		            	FindOccluPair(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, SurfFlag); % detect occlsion (use Ray as detector not the surf Features)
				% data Struction define:
				% PointPix PointDepth FaceSetPickedIND POriReprojM FieldOccluPix OccluDist OccluedFaceSetIND
				% -- allof the size that single ray pass through both FaceSet
				% FaceSetPickedIND1 used in DepthMap size or surfFeature size
 				% OccluedFaceSetIND1 used in DepthMap size
			        if SurfFlag
        				save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '_OccluDetect_Match.mat'],...
	        		        'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1',...
						'OccluedFaceSetIND1','OccluedFaceSetIDRemained1', ...
        	        		'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2',...
						'OccluedFaceSetIND2','OccluedFaceSetIDRemained2');
		        	else    
				        save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '_OccluDetect_Ray.mat'],...
	        		        'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1',...
						'OccluedFaceSetIND1','OccluedFaceSetIDRemained1', ...
        	        		'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2',...
						'OccluedFaceSetIND2','OccluedFaceSetIDRemained2');
                			%'PointPix1','PointDepth1','FaceSetPickedIND1','POriReprojM1','FieldOccluPix1','OccluDist1','OccluedFaceSetIND1',...
                			%'PointPix2','PointDepth2','FaceSetPickedIND2','POriReprojM2','FieldOccluPix2','OccluDist2','OccluedFaceSetIND2');
			        end    
			else
        			if SurfFlag
	        			load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '_OccluDetect_Match.mat']);
			        else    
        				load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '_OccluDetect_Ray.mat']);
	        		end    
		    	end        
	
			% Define occlusion if OccluDist > CentainThreshold
			if defaultPara.Flag.FlagRecipicalOccluDetection
				Mask1 = PointDepth1./OccluDist1 > OcclusionRatioThre | OccluDist1./PointDepth1 > OcclusionRatioThre;
				Mask2 = PointDepth2./OccluDist2 > OcclusionRatioThre | OccluDist2./PointDepth2 > OcclusionRatioThre;
			else
				Mask1 = PointDepth1./OccluDist1 > OcclusionRatioThre;%defaultPara.OccluDistThre;
				Mask2 = PointDepth2./OccluDist2 > OcclusionRatioThre;%defaultPara.OccluDistThre;
			end

			% Action: find new matches or decide plane to remove given the occlusion infomation
			switch lower(Type)
			case 'matches'
				[matches]=MatchPointsGivenOcclusion(defaultPara, ImgInfo1, ImgInfo2, ImgScale1, ImgScale2, Img1, Img2, Img1Index, Img2Index, ...
					Pair,  GlobalScale, Wrlname, PostFixStrAfter, ...
					POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), FaceSetPickedIND1(:,Mask1), ...
					POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), FaceSetPickedIND2(:,Mask2) );	
            case 'corrmatches'
                    if ~isempty(Pair.lamda)
                        depthratio = Pair.lamda(1,:)./Pair.lamda(2,:);
                        % Min add to remove outliers (should be already removed when doing PoseEst.m)
                        Inliers = depthratio > depthratioMin & depthratio < depthratioMax;
                        depthratio = depthratio(Inliers);
                        maxRatio = max(depthratio);
                        minRatio = min(depthratio);
                    else
                        maxRatio = max([GlobalScale(1)/GlobalScale(2) GlobalScale(2)/GlobalScale(1)]);
                        minRatio = 1/maxRatio;
                    end
                    [Matches1 CoeffM1 Inliers1]=CorrolationMatch( defaultPara, Pair, I1, I2, PointPix1(:,Mask1), POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), [minRatio maxRatio]);
		    Pair2_1.R = Pair.R';
	   	    Pair2_1.T = -Pair.R'*Pair.T;
		    if defaultPara.Flag.FlagRefinementDisp
			disp('Start Reverse CorrolationMatches Hold on');
			whos
		    end
                    [Matches2 CoeffM2 Inliers2]=CorrolationMatch( defaultPara, Pair2_1, I2, I1, PointPix2(:,Mask2), POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), [1/maxRatio 1/minRatio]);
		    if defaultPara.Flag.FlagRefinementDisp
			disp('CorrolationMatch.m Finished');
		    end
		    Matches1 = Matches1(:,Inliers1);
		    Matches2 = Matches2(:,Inliers2);
		    CoeffM1 = CoeffM1(:,Inliers1);
		    CoeffM2 = CoeffM2(:,Inliers2);
			 
		    % Check if the Matches are not mutual discard the one with less Coeff(Cross-Corrolation value) ===============
		    Matches = [ Matches1 [Matches2(3:4,:); Matches2(1:2,:)]];
		    CoeffM = [ CoeffM1 CoeffM2];
			% Min used different algorithm than SurFeature Matches
			[Inliers] = CleanMatch(Matches, CoeffM(1,:)); % choose the matches with higher Coeff is the matches is not mutual
			[InliersReverse] = CleanMatch(Matches(:,Inliers), CoeffM(1,Inliers)); % choose the matches with higher Coeff is the matches is not mutual
			Inliers = Inliers(InliersReverse);
			Matches = Matches(:,Inliers);
			CoeffM = CoeffM(:,Inliers);
            
		    if defaultPara.Flag.FlagRefinementDisp
			figure; plotmatches(I1,I2,Matches(1:2,:), Matches(3:4,:),repmat(1:size(Matches,2),2,1), 'Stacking','v','Interactive', 3);
		    end

                    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter 'BeforeFiltering.mat'],'Matches','CoeffM');

  		    % use Coeff as threshould to filter out error matches
		    CoeffMask = CoeffM(1,:) > defaultPara.CoeffMThre;
		    [inlier, Residual] = EpipoPrune(defaultPara, Pair, Matches, ImgScale1);
		    EpipolarResidualMask = Residual < defaultPara.ResidualThre;
		    CoeffRationMask = CoeffM(2,:)./CoeffM(1,:) < defaultPara.coeffratioThre;
		    CoeffRationMask( CoeffM(1,:) == 0) = false;% CoeffM(1,:) == 0 means not satisfy the epipolar constrain
		    Mark = CoeffMask & CoeffRationMask & EpipolarResidualMask;

                    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'Matches','CoeffM','Mark');

		    % =======================================================================
                    if ~isempty(Matches(:,Mark))
                        tempf1 = Matches(1:2,Mark);
                        tempf2 = Matches(3:4,Mark);
                        x_calib = [ inv(defaultPara.InrinsicK1)*[ tempf1; ones(1,size(tempf1,2))];...
                        inv(defaultPara.InrinsicK2)*[ tempf2; ones(1,size(tempf2,2))]];
                        [ lamda1 lamda2 Error] = triangulation( defaultPara, Pair.R, Pair.T, x_calib);
                        % notice lamda re-scale to local model scale
                        lamda1 = lamda1./GlobalScale(1);
                        lamda2 = lamda2./GlobalScale(2);
                        % Storage the match result for later ReInference
                        AddMatch2Model(defaultPara, Wrlname, lamda1, Matches(1:2,Mark), ImgInfo1, ImgScale1, Img1Index, Img2Index, PostFixStrAfter, Error);
                        AddMatch2Model(defaultPara, Wrlname, lamda2, Matches(3:4,Mark), ImgInfo2, ImgScale2, Img2Index, Img1Index, PostFixStrAfter, Error);
			% Storage the New Matches                   
			if defaultPara.Flag.FlagRefinementDisp
				disp('Storaging Occlusion Surf Features Matches');
				save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'Matches','CoeffM','Error','Mark');
			end
                    end

                    
		    case 'remove'
			RemovePlaneGivenOcclusion(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, Img1, Img2, Wrlname, PostFixStrAfter,...
				FaceSetPickedIND1, OccluedFaceSetIND1, FaceSetPickedIND2, OccluedFaceSetIND2, Mask1, Mask2);
		    end
		end
	end
	TimeForOccluDetectionFirstStep = toc;
	if defaultPara.Flag.FlagRefinementDisp
		disp(['Finish OccluDetectionFirstStep in ' num2str(TimeForOccluDetectionFirstStep) ' seconds']);
	end

end % end of if defaultPara.Flag.FlagOccluDetect

% 2) Render of each image individiually -----------------------------------------------------------------------------------
if defaultPara.Flag.FlagRefinementReInference

	tic; % timing the Refinement Step

	for i = 1:NumImg

		switch lower(Type)
		case 'remove'	
			% Data preparing and ReInfernece
    			Default.Wrlname{1} = [Wrlname '_' ImgList{i} '_' PostFixStrAfter];   
	    		Default.filename{1} = ['../' ImgList{i}];
    			load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat']);
    			Default.OutPutFolder = [ defaultPara.OutPutFolder ImgList{i} '/'];
	    		R = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'R');
    			T = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'T');
        
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

		case {'matches','corrmatches'}
			% Data preparing and ReInfernece
    			Default.Wrlname{1} = [Wrlname '_' ImgList{i} '_' PostFixStrAfter];   
	    		load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat']);
    			R = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'R');
    			T = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'T');
	    		Scale = LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'Scale');
			SingleImgReInference(defaultPara, model, ImgList{i}, LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'GroundLevel'),Wrlname, ...
				R, T, Scale ,PostFixStrAfter);
		end

		%   Build Meta Wrl file =========
	  	Path = [ defaultPara.OutPutFolder Wrlname '_' PostFixStrAfter '.wrl'];
   		InLinePath = ['./'  ImgList{i} '/' Default.Wrlname{1} '.wrl'];
	    	Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
    		if any(isnan(Q))
        		Q = zeros(4,1);
	    	end
    		WRLT = T;
	    	WRLT(3) = -WRLT(3);% for Wrl (-) z component;
    		Scale = LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'Scale');
	    	if i == 1
        		BuildVrmlMetaModel(1, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1));
	    	else
        		BuildVrmlMetaModel(0, defaultPara.OutPutFolder, Path, InLinePath, Q, WRLT, repmat(Scale,3,1)); 
	    	end
    		% ===============================
	
	end
	
	TimeForRefinementStep = toc;
	if defaultPara.Flag.FlagRefinementDisp
		disp(['Finish RefinementStep in ' num2str(TimeForRefinementStep) ' seconds']);
	end
end % end of if defaultPara.Flag.FlagRefinementReInference

return;
