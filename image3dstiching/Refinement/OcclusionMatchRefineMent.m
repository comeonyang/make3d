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
function [] = OcclusionMatchRefineMent(defaultPara, Wrlname, PairList)

% This function do three thing
% 1) detect the occlusion from the current model
% 2) Match the Occlued area to find any match 
% 3) reinference including the new matches information

% Input:
% 	1) defaultPara - all parameter share for the whole stitching3d folder
%	2) Wrlname - the name of the model, different by the PairList
%	3) PairList - the linearly order of image that been added to the model

% Fixed constant --- need to be put in defaultPara by Min later
NegI = diag([1 1 -1]);
PostFixStr = 'NonMono'
PostFixStrAfter = 'NonMonoOccluMatched'
% 1) Run occlusion detection and Matcheing for every pair of images

ImgList = unique(PairList);
NumImg = length(ImgList);

for i = 1:NumImg
	for j = (i+1):NumImg

%  	detect occlusion for a pair of image

	[ImgInfo1 ImgInfo2 Pair GlobalScale] = ... % What to do with Pairs haven't been matched ???????????????????
		LoadDataForFindOcclu(defaultPara, Wrlname, ImgList{i}, ImgList{j}, ...
                		     0, PostFixStr); % load all imformation needed for occlusion detection

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
        [Region1 PointPix1 POriReprojM1 PoccluM1 OccluDist1 Region2 PointPix2 POriReprojM2 PoccluM2 OccluDist2] = ...
                    FindOccluPair(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, 1); % detect occlsion
        save([ defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect.mat'],...
                'Region1','PointPix1','POriReprojM1','PoccluM1','OccluDist1','Region2','PointPix2','POriReprojM2','PoccluM2','OccluDist2');
    else
        load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_OccluDetect.mat']);
    end        

	% Define occlusion if OccluDist > CentainThreshold
	Mask1 = OccluDist1 > defaultPara.OccluDistThre;
	Mask2 = OccluDist2 > defaultPara.OccluDistThre;
	[Region1 PointPix1 POriReprojM1 PoccluM1 Region2 PointPix2 POriReprojM2 PoccluM2] = ...
		OccluPruning(Mask1, Mask2, Region1, PointPix1, POriReprojM1, PoccluM1, OccluDist1, ...
			     Region2, PointPix2, POriReprojM2, PoccluM2, OccluDist2);

	
%  	match occlusion part for a pair of image 
	% Prepare the Constrain to run the matching 
	[ Rc1 ConS1 ConSRough1] = EndPoint2BoxConS(defaultPara, H, V, POriReprojM1, PoccluM1);
        [ Rc2 ConS2 ConSRough2] = EndPoint2BoxConS(defaultPara, H, V, POriReprojM2, PoccluM2);
	Vector2Ipoint([Rc1; ConS1],[defaultPara.Fdir '/surf/'],['RConS_' Img1]);
        Vector2Ipoint([Rc2; ConS2],[defaultPara.Fdir '/surf/'],['RConS_' Img2]);
        Vector2Ipoint([ConSRough1],[defaultPara.Fdir '/surf/'],['RConSRough_' Img1]);
        Vector2Ipoint([ConSRough2],[defaultPara.Fdir '/surf/'],['RConSRough_' Img2]);
 	tic;
        cd match
%                 system(['./surfOccluMatch.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' OccluDense ' '0.1 0.2']);    % Parameter still need to be changed//Min
        cd ..
        toc
        [f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, [ defaultPara.Type 'OccluDense'], 1, 1);
        figure(200); plotmatches(I1,I2,f1, f2,matches, 'Stacking','v','Interactive', 3);
	saveas(200,[ defaultPara.ScratchFolder Img1 '_' Img2 '_OccluMatches'],'jpg');

% 	Process the matches
	% Pruning by epipolarline
    [inlier] = EpipoPrune(defaultPara, Pair, [f1(:,matches(1,:)); f2(:,matches(2,:))], (ImgScale1+ImgScale2)/2);
    matches = matches(:,inlier);
    figure(201); plotmatches(I1,I2,f1, f2,matches, 'Stacking','v','Interactive', 3);
	saveas(201,[ defaultPara.ScratchFolder Img1 '_' Img2 '_OccluMatchesPrune'],'jpg');

    
	% Triangulation
    if ~isempty(matches)
    	tempf1 = f1(:,matches(1,:));
        tempf2 = f2(:,matches(2,:));
    	x_calib = [ inv(defaultPara.InrinsicK1)*[ tempf1; ones(1,size(tempf1,2))];...
                 inv(defaultPara.InrinsicK2)*[ tempf2; ones(1,size(tempf2,2))]];
    	[ lamda1 lamda2 Error] = triangulation( defaultPara, Pair.R, Pair.T, x_calib);
        % notice lamda re-scale to local model scale
        lamda1 = lamda1./GlobalScale(1);
        lamda2 = lamda2./GlobalScale(2);
    	% Storage the match result for later ReInference
        AddMatch2Model(defaultPara, Wrlname, lamda1, f1(:,matches(1,:)), ImgInfo1, ImgScale1, i, j, PostFixStrAfter);
    	AddMatch2Model(defaultPara, Wrlname, lamda2, f2(:,matches(2,:)), ImgInfo2, ImgScale2, i, j, PostFixStrAfter);
    end    

%  	storage the new Triangulated information
	end
end

% 2) ReInference of each image individiually

for i = 1:NumImg
	
%	Data preparing and ReInfernece
    Default.Wrlname{1} = [Wrlname '_' ImgList{i} '_OccluMatches'];   
    load([defaultPara.Fdir '/data/' ImgList{i} '/' Wrlname '_' ImgList{i} '_' PostFixStrAfter '.mat']);
    R = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'R');
    T = LoadModelStatus( defaultPara.Fdir, Wrlname, ImgList{i}, 'T');
    Scale = LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'Scale');
	SingleImgReInference(defaultPara, model, ImgList{i}, LoadModelStatus(defaultPara.Fdir, Wrlname, ImgList{i}, 'GroundLevel'),Wrlname, ...
	R, T, Scale ,PostFixStrAfter);
    
%   Build Meta Wrl file

    Path = [ defaultPara.OutPutFolder Wrlname 'Occlu.wrl'];
    InLinePath = ['./'  ImgList{i} '/' Default.Wrlname{1} '.wrl'];
    % Possibly to be wrong ==========
    %       Angle = recoverAlphasFromU(reshape(R,1,[]));
    %       Q = GetQauternionFrom2Rotation(zeros(3,1), Angle, false);%[1 0 0 0]';
    Q = Rotation2Q(NegI*R*NegI);% for Wrl (-) z component;
    if any(isnan(Q))
        Q = zeros(4,1);
    end
    % ===============================
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
