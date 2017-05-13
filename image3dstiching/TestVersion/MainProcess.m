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
function []= MainProcess(defaultPara, Fdir);

% This function is the main function tha handle:
% 1) Extracting info mation : Exif GPS IMU
% 2) MonoDepth estimation
% 3) Ransac
% 4) MetricReconstruction
% 5) Re-inference
% 6) Rendering

% initialize parameters
NotAllSurfCompute = 0;
Type = '_RConS';
FlagInitialPathLearningCode = 0;
FlagRenderMonoWrlJintly = 1;
ops = sdpsettings('solver','sedumi','verbose',1);
% parameter for Mono depth Estimation ------------------------
ScratchFolder = '/afs/cs/group/reconstruction3d/scratch/temp';
OutPutFolder = '/afs/cs/group/reconstruction3d/scratch/3DmodelMultipleImage/';
ParaFolder = '/afs/cs/group/reconstruction3d/scratch/Para/';
taskName = '';
Flag.DisplayFlag = 0;
Flag.IntermediateStorage = 0;
Flag.FeaturesOnly = 0;
Flag.NormalizeFlag = 1;
Flag.BeforeInferenceStorage = 0;
Flag.NonInference = 0;
Flag.AfterInferenceStorage = 1;
% ------------------------------------------------------------

% initialize variables
ImgInfo = [];
ImgAddedList = [];

% 1) Extracting info
if NotAllSurfCompute 
	system(['./surf/surfFeatureDir.sh ' Fdir]);
end

[ImgInfo] = ExifExtractDir(Fdir,ImgInfo); % Exif info
[ImgInfo] = ExtractGPSInfo(defaultPara, Fdir, ImgInfo); % GPS info
[ImgInfo] = ExtractRotationInfo( defaultPara, Fdir, ImgInfo); % IMU info

% start loop by getting input
while length(ImgAddedList) <= length(ImgInfo)
	
	% request user input
	NewImg = input('add image', 's');

	% 2) MonoDepth estimation for New Imag------------------
	NewImgPath = [Fdir '/jpg/' NewImg '.jpg'];
	if system(['ls ' ScratchFolder '/' NewImg '__AInfnew.mat']);
		cd ../LearningCode
		if ~FlagInitialPathLearningCode
			InitialPath;
			FlagInitialPathLearningCode = 1;
 	 	end
		OneShot3dEfficient(NewImgPath, OutPutFolder,...
	        taskName,...% taskname will append to the imagename and form the outputname
        	ScratchFolder,... % ScratchFolder
	        ParaFolder,...
        	Flag...  % All Flags 1) intermediate storage flag
        	);
      		cd ../multipleImages
   	end
   	load( [ScratchFolder '/' NewImg '__AInfnew.mat']);

        % rename and clear vriables
	Depth2 = FitDepth;
        Sup2 = Sup;
	clear MedSup MultiScaleSupTable Ray RayOri SupNeighborTable SupOri depthMap maskG maskSky Sup FitDepth;
	% ------------------------------------------------------

	% skip the following step if ImgAddedList is empty
	if isempty(ImgAddedList)
		Depth = Depth2;
		Sup = Sup2;
		save([Fdir '/data/' NewImg '_Data.mat'],'Depth','Sup');
	        ImgAddedList = [ImgAddedList NewImg];
		continue;
	end

	% Use track concept to define track
	PairImg = input('Pair image in addedList', 's');

	% build up data need for processing
	BuildDataPairImg;

	% 3) Ransac --------------------------------------------
	% Estimate Rotation and translation
	[R T] = InitPoseMeas(defaultPara, ImgInfo1, ImgInfo2);
	T([2 5]) = 0;
    
	% Generate Constrain for SurfMatch searhing region
	[ Rc1, Rc2, ConS1, ConS2]=CalMatchSearchRegin(defaultPara, R, T, I1, I2, f1, f2, D1, D2, 1, 1);
	Vector2Ipoint([Rc1; ConS1],[Fdir '/surf/'],['RConS_' PairImg]);	
	Vector2Ipoint([Rc2; ConS2],[Fdir '/surf/'],['RConS_' NewImg]);	

	% Constrasin Search Surf Match
%	if system(['ls ' Fdir '/surf_matches/' PairImg '-' NewImg '.match_RConS']) ...
%	   && system(['ls ' Fdir '/surf_matches/' NewImg '-' PairImg '.match_RConS'])
		cd match
		system(['./surfMatchRConS.sh ' Fdir ' ' PairImg ' ' NewImg]);    
		cd ..
%	end
	% general Ransac
	[f1, f2, matches] = readSurfMatches(PairImg, NewImg, Fdir, Type, 0, 1);
% 	displaySurfMatches(Fdir, PairImg, NewImg, Type, 0, 1);
    
	[DV DH] = size(Depth2);
	[IV IH] = size(I2);
	Scale = [[DV DH];[IV IH]];
	[IND2]=ReScaleInd2Sub(Scale,f2(:,matches(2,:)));
	D2 = Depth2(IND2);
	[DV DH] = size(Depth1);
	[IV IH] = size(I1);
	Scale = [[DV DH];[IV IH]];
	[IND1]=ReScaleInd2Sub(Scale,f1(:,matches(1,:)));
	D1 = Depth1(IND1);
	[F, inliers, NewDist, fail]=GeneralRansac(defaultPara, f1, f2, matches, D1, D2);
   	figure;
    	plotmatches(I1,I2,f1, f2,matches(:,inliers), 'Stacking', 'v', 'Interactive', 2);

    % maually remove outliers and re-estimate F    
    matches = matches(:,inliers);
    [F, inliers, fail] = ransacmatches(defaultPara, f1, f2, matches, I1, I2, 1);
    
	% Nonlinear Refinement
%  	[F, newinliers, fail] = RansacNonlinearReFinement(f1, f2, matches, F, inliers, I1, I2, 1);
	% ------------------------------------------------------

	% 4) MetricReconstruction ------------------------------
    % Initial depth estimated from fundamatal matrix
    E = defaultPara.InrinsicK2*F*defaultPara.InrinsicK1;
    x = [ inv(defaultPara.InrinsicK1)*[ f1(:,matches(1,inliers)); ones(1,length(inliers))];...
          inv(defaultPara.InrinsicK2)*[ f2(:,matches(2,inliers)); ones(1,length(inliers))]];
    [ R_f, T_f, lamda1, lamda2, inlier] = EstPose( E, x, R(1:3,:), T(1:3), false);
    
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.  Note that 'fundmatrix' will also call
    % 'normalise2dpts' but the code in 'ransac' that calls the distance
    % function will not - so it is best that we normalise beforehand.
    x1 = defaultPara.InrinsicK1*x(1:3,inlier);
    x2 = defaultPara.InrinsicK2*x(4:6,inlier);
    [x1_normalized, T1] = normalise2dpts(x1);
    [x2_normalized, T2] = normalise2dpts(x2);
    
    % Projective resonstruction
    [ P, X, Outliner] = ProjectionFactorization( 1, 1, ...
                        cat(3, x1_normalized, x2_normalized), cat(3,lamda1(inlier)', lamda2(inlier)'));
	%    [ P, X, Outliner] = ProjectionFactorization_missing_data...
	%                ( ARes, BRes, x_im, inc, depth);

    % Find out the H (projective transformation)
    matches = matches(:,inliers);
    x1 = [f1(:,matches(1,inlier)); ones(1,length(inlier))];
    x2 = [f2(:,matches(1,inlier)); ones(1,length(inlier))];
    Scale = [[DV DH];[IV IH]];
    [IND1]=ReScaleInd2Sub(Scale,f1(:,matches(1,inlier)));
    D1 = Depth1(IND1);
    Scale = [[DV DH];[IV IH]];
    [IND2]=ReScaleInd2Sub(Scale,f2(:,matches(1,inlier)));
    D2 = Depth2(IND2);
    [R_mr T_mr D1_modified D2_modified] = MetricReconEstPose(P, R(1:3,:), T(1:3), T1, T2, D1, D2, x1, x2, defaultPara);       
    
    % Triangulate the depth
if false
    x1 = inv(T1)*x1;
    x2 = inv(T1)*x2;
    x1(:,Outliner) = [];     
    x2(:,Outliner) = [];     
    X = [x1; x2];
    [lamda1, lamda2] = Triangulate( R_mr, T_mr, X);
    % Rescale the Depth1 and Depth2
    x1 = defaultPara.InrinsicK1*x1; 
    x2 = defaultPara.InrinsicK2*x2;
    [DV DH] = size(Depth1);
    [IV IH] = size(I1);
    Scale = [[DV DH];[IV IH]];
    [IND1]=ReScaleInd2Sub(Scale,x1);
    D1 = Depth1(IND1);
    [DV DH] = size(Depth2);
    [IV IH] = size(I2);
    Scale = [[DV DH];[IV IH]];
    [IND2]=ReScaleInd2Sub(Scale,x2);
    D2 = Depth2(IND2);
    Scale1 = sdpvar(1,1);
    Scale2 = sdpvar(1,1);
    F = set(Scale1>=0)+set(Scale2 >=0);
    sol = solvesdp(F,norm(lamda1 - Scale1*D1', 1)+ norm(lamda2 - Scale2*D2', 1),ops);
    Scale1 = double(Scale1); 
    Scale2 = double(Scale2); 
    Depth1 = Depth1*Scale1;
    Depth2 = Depth2*Scale2;
else
	Scale1 = median(D1_modified./D1);
	Scale2 = median(D2_modified./D2);
        Depth1 = Depth1*Scale1;
    	Depth2 = Depth2*Scale2;
%	Depth1(IND1) = D1_modified;
%	Depth2(IND2) = D2_modified;
end	
    
	% ------------------------------------------------------	
   
	% Test R_mr and T_mr to show the Mono-VRml jointly  according to R_mr and T_mr
	if FlagRenderMonoWrlJintly
		[WrlPosition1 PositionTex1] = PreWrlData(defaultPara.InrinsicK1, Depth1, I1, eye(3), zeros(3,1));
		WrlFacestHroiReduce(WrlPosition1, PositionTex1, Sup1, PairImg, [PairImg  '_' NewImg '_MonoJointly'], OutPutFolder, 0, 0);
		system(['cp ' Fdir '/jpg/' PairImg '.jpg '  OutPutFolder PairImg '.jpg']);
		[WrlPosition2 PositionTex2] = PreWrlData(defaultPara.InrinsicK2, Depth2, I2, R_mr, T_mr);
		WrlFacestHroiReduce(WrlPosition2, PositionTex2, Sup2, NewImg, [PairImg  '_' NewImg '_MonoJointly'], OutPutFolder, 0, 1);
		system(['cp ' Fdir '/jpg/' NewImg '.jpg ' OutPutFolder NewImg '.jpg']);
	end
 
    % 5) Reinfererce
    
	ImgAddedList = [ImgAddedList NewImg];

    
	% save Data
	Depth = Depth2;
	Sup = Sup2;
	save([Fdir '/data/' NewImg '_Data.mat'],'Depth','Sup');
end
