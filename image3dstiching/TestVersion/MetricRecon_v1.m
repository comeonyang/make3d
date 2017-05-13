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
function [] = MetricRecon(defaultPara, ImgInfo, Fdir, left, right, Wrlname, appendOpt, RenderFlag)
% This function estimate the orientation and translation of a pair of imaegs
% Always set the left image as the world coordinate
% assumption:
% 1) image center is in the cneter of the pixel coordinate (e. round(2272/2) round(1704/2))
% 2) No Skew
% 3) Aspect ratio is know
% processure:
% 1) run the ransac matching 
% 2) retrive the camera matrix(P) and structure(X) up to projective transformation
% 3) Using the assumption of the camera intrinsic parameter to do linearized self-calibration
% 4) rescale the learned depth of two images
% 5) plane parameter inference of a pair of images

% Setup parameters
Type = '_RConS';
ResultFolder = '/afs/cs/group/reconstruction3d/scratch/3DmodelMultipleImage';
EstimateIntrinsicCameraPara % setup the intrinsic camera parameters
image3dstichingSetupUp% setup for the new image3d stitching functinos
ACalibratedFlag = 1; % set to 1 if using EstimateIntrinsicCameraPara
BCalibratedFlag = 1; % set to 1 if using EstimateIntrinsicCameraPara
UnDoAspecRatioFlag = 0;
ClosestMatchfindFlag = 0;
RenderWrlDirectlyFlag = 1;
GroundStaticFlag = 0;
GroundStaticWeight = 5;
solverVerboseLevel = 0;
Match3DThreshould = 5; % meters
opt = sdpsettings('solve','sedumi','cachesolvers',1,'verbose',solverVerboseLevel);

AImgPath = [Fdir '/jpg/' left '.jpg'];
BImgPath = [Fdir '/jpg/' right '.jpg'];
OutPutFolder = '/afs/cs/group/reconstruction3d/scratch/3DmodelMultipleImage/'; 
taskName = '';%(Not Used) taskname will append to the imagename and form the outputname
Flag.DisplayFlag = 0;
Flag.IntermediateStorage = 0;
Flag.FeaturesOnly = 0;
Flag.NormalizeFlag = 1;
Flag.BeforeInferenceStorage = 0;
Flag.NonInference = 0;
Flag.AfterInferenceStorage = 1;
%ScratchFolder = '/afs/cs/group/reconstruction3d/scratch/IMStorage'; % ScratchFolder
ScratchFolder = '/afs/cs/group/reconstruction3d/scratch/temp'; % ScratchFolder
ParaFolder = '/afs/cs/group/reconstruction3d/scratch/Para/'; % default parameter folder


% load images
imgA = imread([ Fdir '/jpg/' left '.jpg']);
AimgCameraParameters = exifread( [ Fdir '/jpg/' left '.jpg']);
[Ay Ax dummy] = size(imgA);
imgB = imread([ Fdir '/jpg/' right '.jpg']);
BimgCameraParameters = exifread( [ Fdir '/jpg/' left '.jpg']);
[By Bx dummy] = size(imgB);

% copy image to the Outputfolder
system(['cp ' Fdir '/jpg/' left '.jpg ' OutPutFolder]);
system(['cp ' Fdir '/jpg/' right '.jpg ' OutPutFolder]);

% known assumption
% 1. camera center
if ACalibratedFlag
   AUx = Default.Ox; % in pixel coordinate % Default round(Px/2)
   AUy = Default.Oy; % in pixel coordinate % Default round(Py/2)
elseif 0
   AUx = 1;
   AUy = 1;
else
   AUx = round(Ax/2); % in pixel coordinate % Default round(Px/2)
   AUy = round(Ay/2); % in pixel coordinate % Default round(Py/2)
end
if BCalibratedFlag
   BUx = Default.Ox; % in pixel coordinate % Default round(Px/2)
   BUy = Default.Oy; % in pixel coordinate % Default round(Py/2)
elseif 0
   BUx =1;
   BUy =1; 
else
   BUx = round(Bx/2); % in pixel coordinate % Default round(Px/2)
   BUy = round(By/2); % in pixel coordinate % Default round(Py/2)
end

% 2. camera aspect ratio
if ACalibratedFlag
   AAspectRatio = Default.fx/Default.fy; % in pixel coordinate % Default round(Px/2)
elseif 0
   AAspectRatio = 1;
else
   AAspectRatio = 1; % in pixel coordinate % Default round(Px/2)
end
if BCalibratedFlag
   BAspectRatio = Default.fx/Default.fy; % in pixel coordinate % Default round(Px/2)
elseif 0
   BAspectRatio = 1;
else
   BAspectRatio = 1; % in pixel coordinate % Default round(Px/2)
end

% 1) run the matching to get F and match point
if system(['ls ' Fdir '/surf/*.surf'])
  cd ./match/
%  system(['./surfFeaturesAndMatchesDir.sh ' Fdir]);
  system(['./surfMatch.sh ' Fdir '' left '' right]);
  cd ..
end

I1=imreadbw([Fdir '/pgm/' left '.pgm']);
I2=imreadbw([Fdir '/pgm/' right '.pgm']);


[f1, f2, matches] = readSurfMatches(left, right, Fdir, Type, false, true);
% displaySurfMatches(Fdir, left, right, Type, 0, 1);
% ============Load mono info==========
% process imageA
% cd ../LearningCode
% InitialPath
% cd ../multipleImages
if appendOpt
   load( [ScratchFolder '/' left '_NonMono.mat']);
else
   if system(['ls ' ScratchFolder '/' left '__AInfnew.mat']);
      cd ../LearningCode
      OneShot3dEfficient(AImgPath, OutPutFolder,...
        taskName,...% taskname will append to the imagename and form the outputname
        ScratchFolder,... % ScratchFolder
        ParaFolder,...
        Flag...  % All Flags 1) intermediate storage flag
        );
      cd ../multipleImages
   end
   load( [ScratchFolder '/' left '__AInfnew.mat']);
   % ==================
   ARotation = eye(3);
   ATranslation = zeros(3,1);
   AHistory=[];
   %AdepthMapRaw = depthMap;
   AdepthMap = full(FitDepth);
   ARayAll = Ray;
   ARayOri = RayOri;
   ASup = Sup;
   AMedSup = MedSup;
   ASupOri = SupOri;
   AMedSup = MedSup;
   ASupNeighborTable =SupNeighborTable;
   AmaskSky = maskSky;
   AmaskG = maskG;
   AMultiScaleSupTable = MultiScaleSupTable;
   Aconstrain.RayMatched = [  ];
   Aconstrain.Depth_modified = [  ];
   Aconstrain.SupMatched = [ ];
% =================
end
% process imageB
if system(['ls ' ScratchFolder '/' right '__AInfnew.mat']);
   cd ../LearningCode
   OneShot3dEfficient(BImgPath, OutPutFolder,...
     taskName,...% taskname will append to the imagename and form the outputname
     ScratchFolder,... % ScratchFolder
     ParaFolder,...
     Flag...  % All Flags 1) intermediate storage flag
     );
   cd ../multipleImages
end
load([ScratchFolder '/' right '__AInfnew.mat']);
% ==================
%BdepthMapRaw = depthMap;
BHistory=[];
BdepthMap = full(FitDepth);
BRayAll = Ray;
BRayOri = RayOri;
BSup = Sup;
BMedSup = MedSup;
BSupOri = SupOri;
BMedSup = MedSup;
BSupNeighborTable =SupNeighborTable;
BmaskSky = maskSky;
BmaskG = maskG;
BMultiScaleSupTable = MultiScaleSupTable;
% =================
% robut estimate of matches
IND1 = sub2ind([55 305],max(min(round(f1(2,matches(1,:))/size(I1,1)*55),55),1),max(min(round(f1(1,matches(1,:))/size(I1,2)*305),305),1));
IND2 = sub2ind([55 305],max(min(round(f2(2,matches(2,:))/size(I1,1)*55),55),1),max(min(round(f2(1,matches(2,:))/size(I1,2)*305),305),1));
D1 = AdepthMap(IND1);
D2 = BdepthMap(IND2);
[F, newinliers, NewDist, fail]=GeneralRansac(defaultPara, f1, f2, matches, D1, D2);
figure;
    	plotmatches(I1,I2,f1, f2,matches(:,newinliers), 'Stacking', 'v', 'Interactive', 2);
%  [ F, newinliers, fail]=RansacOnPairMatches(defaultPara, f1, f2, matches, I1, I2, D1, D2, 1)
% =================
% ====================================

AMatch = [f1(1,matches(1,newinliers)) ; f1(2,matches(1,newinliers))]; % [x ; y] pixel coordinate (left top as origin)
BMatch = [f2(1,matches(2,newinliers)) ; f2(2,matches(2,newinliers))]; % [x ; y] pixel coordinate (left top as origin)
% BMatch = [f1(1,matches(1,newinliers)) ; f1(2,matches(1,newinliers))]; % [x ; y] pixel coordinate (left top as origin)
% AMatch = [f2(1,matches(2,newinliers)) ; f2(2,matches(2,newinliers))]; % [x ; y] pixel coordinate (left top as origin)
% =================================================================

% ==========undo the camera center and the camera aspect ratio if needed ====================
% to improve the condition of the camera matrix
% translate the coordinate into camera center coordinate
AMatchCenter = [ [AMatch(1,:) - AUx];...
                 [Ay + 1 - AMatch(2,:) - AUy] ]; % [x ; y] pixel coordinate (center as origin)
BMatchCenter = [ [BMatch(1,:) - BUx];...
                 [By + 1 - BMatch(2,:) - BUy] ]; % [x ; y] pixel coordinate (center as origin)
if UnDoAspecRatioFlag 
   AScale = [ Ax; Ay]; % important
   BScale = [ Bx; By]; % important
else
   AScale = [ Ax; Ax]; % important
   BScale = [ Bx; Bx]; % important
end
% scaling to improve condition
AMatchCenter_normalized = AMatchCenter./repmat(AScale,1,size( AMatchCenter,2));
BMatchCenter_normalized = BMatchCenter./repmat(BScale,1,size( BMatchCenter,2));
% =================================================================
Ray = permute(Ray, [ 2 3 1]);
[APositionAll] = im_cr2w_cr( AdepthMap,Ray);
[Dy Dx] = size(AdepthMap);
AMatchDepthRes = [ ( AMatch(1,:)-0.5)/Ax*Dx+0.5 ;...
                   ( AMatch(2,:) -0.5)/Ay*Dy+0.5]; % [x ; y] in Dy Dx Resolution
% the closet approximate of the match points to the depthMap index
Aindex = min(max( round( AMatchDepthRes(1,:))-1, 0), Dx)*Dy+ min(max( round( AMatchDepthRes(2,:)), 0), Dy);
ADepthMatch = AdepthMap( Aindex);
ASampleImCoordYSmall = (( Dy+1-AMatchDepthRes(2,:))-0.5)/Dy - Default.Oy_default; 
ASampleImCoordXSmall = (AMatchDepthRes(1,:) -0.5)/Dx - Default.Ox_default;
ARayMatch = permute( RayImPosition( ASampleImCoordYSmall, ASampleImCoordXSmall,...
                          Default.a_default, Default.b_default, ...
                          Default.Ox_default,Default.Oy_default), [3 2 1]); %[ 3 horiXSizeLowREs VertYSizeLowREs]
APositionMatch = [ARayMatch'.*repmat( ADepthMatch', 1, 3)];

[BPositionAll] = im_cr2w_cr( BdepthMap,Ray);
BMatchDepthRes = [ ( BMatch(1,:)-0.5)/Bx*Dx+0.5 ; ( BMatch(2,:)-0.5)/By*Dy+0.5]; % [x ; y]
Bindex = min(max( round( BMatchDepthRes(1,:))-1, 0), Dx)*Dy+ min(max( round( BMatchDepthRes(2,:)), 0), Dy);
BDepthMatch = BdepthMap( Bindex);
BSampleImCoordYSmall = (( Dy+1-BMatchDepthRes(2,:))-0.5)/Dy - Default.Oy_default;
BSampleImCoordXSmall = ( BMatchDepthRes(1,:) -0.5)/Dx - Default.Ox_default;
BRayMatch = permute( RayImPosition( BSampleImCoordYSmall, BSampleImCoordXSmall,...
                          Default.a_default, Default.b_default, ...
                          Default.Ox_default,Default.Oy_default), [3 2 1]); %[ 3 horiXSizeLowREs VertYSizeLowREs]
BPositionMatch = (BRayMatch').*repmat( BDepthMatch', 1, 3);

% 2) retrive the camera matrix(P) and structure(X) up to projective transformation
[ P, X, OutlinerCameraCoord] = ProjectionFactorization( [Ax; Ay]./AScale, [Bx; By]./BScale, ...
                                   cat( 3, [ AMatchCenter_normalized(1,:); AMatchCenter_normalized(2,:); ones(1, size(AMatchCenter_normalized ,2)) ],...
                                   [ BMatchCenter(1,:); BMatchCenter(2,:);  ones(1, size(AMatchCenter_normalized ,2))]),...
                                   cat( 3, APositionMatch(:,3)', BPositionMatch(:,3)'));

% remove outliners
Aindex(OutlinerCameraCoord) = []; 
ADepthMatch(OutlinerCameraCoord) = [];
ARayMatch(:,OutlinerCameraCoord) = [];
APositionMatch(OutlinerCameraCoord,:) = [];
ASampleImCoordYSmall(OutlinerCameraCoord) = [];
ASampleImCoordXSmall(OutlinerCameraCoord) = [];
Bindex(OutlinerCameraCoord) = []; 
BDepthMatch(OutlinerCameraCoord) = [];
BRayMatch(:,OutlinerCameraCoord) = [];
BPositionMatch(OutlinerCameraCoord,:) = [];
BSampleImCoordYSmall(OutlinerCameraCoord) = [];
BSampleImCoordXSmall(OutlinerCameraCoord) = [];

% 3) Using the assumption of the camera intrinsic parameter to do linearized self-calibration
[ U S V] = svd(P(1:3,:));
% ============ simple try ================================
 P1 = P(1:3,:);
 P2 = P(4:6,:);
 Ksimple = diag( [ Default.fx/AScale(1) Default.fy/AScale(2) 1]);
 H_a_simple = P1\[Ksimple]; % the upper 4x3 part of the projective transformation matrix
 [ U S V] = svd(P1);
 H_b_simple = V(:,4); % the last 4x1 column vector
 K2simple_square = sdpvar(3,1);
 F = set(diag(K2simple_square) >=0);
 sol = solvesdp(F , norm( P2*H_a_simple*H_a_simple'*P2' - diag(K2simple_square),'fro'), opt);
 K2simple_square = double(K2simple_square);
 K2simple = sqrt(K2simple_square);
 R2_est = diag( 1./K2simple)*P2*H_a_simple; 
if GroundStaticFlag
  disp('Ground Static');
  R2_ground_constrain = sdpvar(3,3);
  sol =solvesdp( [], norm( R2_ground_constrain - R2_est,'fro') + GroundStaticWeight*norm([0 1 0]' - R2_ground_constrain*[0 1 0]'), opt);
  R2_ground_constrain = double( R2_ground_constrain);
  disp('[0; 1; 0] - R2_ground_constrain*[0; 1; 0]');
  norm([0 1 0]' - R2_ground_constrain*[0 1 0]')
  R2 = R2_ground_constrain * (R2_ground_constrain'*R2_ground_constrain)^(-.5);
else
 R2 = R2_est*(R2_est'*R2_est)^(-0.5);
end
 % check if norm([0 1 0]' - R2*[0 1 0]');
 disp('[0; 1; 0] - R2*[0; 1; 0]');
 norm([0 1 0]' - R2*[0 1 0]')
%  pause;

 % set image B's groud and vertical vector;
% BRotation = R2'*ARotation';
 BRotation = ARotation*R2';
 %BDirectionFromReference = R2*eye(3);
 ADirectionFromReference = ARotation'*eye(3);
 BDirectionFromReference = BRotation'*eye(3);

 T_simple_before_R = R2'*inv(diag(sqrt(K2simple_square)))*P2*H_b_simple;
 yalmip('clear');
 % after calibration
 A.a_default = Default.TrainHoriXSize/Default.fx;
 A.b_default = Default.TrainVerYSize/Default.fy;
 B.a_default = A.a_default;%1704/(K2simple(1)*AScale(1)/K2simple(3));
 B.b_default = A.b_default;%2272/(K2simple(2)*AScale(2)/K2simple(3));
 ARay_after_calibration = permute( RayImPosition( ASampleImCoordYSmall, ASampleImCoordXSmall,...
                          A.a_default, A.b_default), [3 2 1]); %[ 3 horiXSizeLowREs VertYSizeLowREs]
 BRay_after_calibration = permute( RayImPosition( BSampleImCoordYSmall, BSampleImCoordXSmall,...
                          B.a_default, B.b_default), [3 2 1]); %[ 3 horiXSizeLowREs VertYSizeLowREs]
 % solve for the translation scaling factor
 T = sdpvar(3,1);
 ADepth_modified = sdpvar(1,length(ADepthMatch));
 BDepth_modified = sdpvar(1,length(BDepthMatch));
 BDepthScale = sdpvar(1);
 sol =solvesdp( [], norm( reshape( (ARay_after_calibration.*repmat( ADepth_modified, 3,1))' - ...
                    ( R2'*( BRay_after_calibration.*repmat( BDepth_modified, 3, 1) ) + repmat(T,1,length(ADepthMatch)) )' ,1,[]),1)+...
                    norm(ADepth_modified - ADepthMatch,1) + norm(BDepth_modified - BDepthScale*BDepthMatch,1), opt);
 T_simple_before_R = double(T);
 ADepth_modified = double(ADepth_modified);
 BDepth_modified = double(BDepth_modified);    
 BDepthScale = double(BDepthScale);
 ADepthScale = 1;

 % clean the 3D matches if the |AMatchPosition - BMatchPosition | > Match3DThreshould =================
 DistMatchError = ARay_after_calibration.*repmat( ADepth_modified, 3,1) - ...
                    ( R2'*( BRay_after_calibration.*repmat( BDepth_modified, 3, 1) ) + repmat(T_simple_before_R,1,length(ADepthMatch)) );
 DistMatchError = norms(DistMatchError);
 OutlinerMark = DistMatchError > Match3DThreshould;

 Aindex(OutlinerMark) = [];
 ADepthMatch(OutlinerMark) = [];
 ARayMatch(:,OutlinerMark) = [];
 APositionMatch(OutlinerMark,:) = [];
 ASampleImCoordYSmall(OutlinerMark) = [];
 ASampleImCoordXSmall(OutlinerMark) = [];
 Bindex(OutlinerMark) = [];
 BDepthMatch(OutlinerMark) = [];
 BRayMatch(:,OutlinerMark) = [];
 BPositionMatch(OutlinerMark,:) = [];
 BSampleImCoordYSmall(OutlinerMark) = [];
 BSampleImCoordXSmall(OutlinerMark) = [];
 ARay_after_calibration(:, OutlinerMark) = [];
 BRay_after_calibration(:, OutlinerMark) = [];
 ADepth_modified( OutlinerMark) = [];
 BDepth_modified( OutlinerMark) = [];

 % ===================================================================================================
 if ClosestMatchfindFlag
    % solve the opt depth for matches of image A and image B
    opt = sdpsettings('solver','sedumi');
    AClosestDepth = sdpvar( 1,length(ADepthMatch));
    BClosestDepth = sdpvar( 1,length(BDepthMatch));
    F = set(AClosestDepth >=0) + set(BClosestDepth >=0);
    sol = solvesdp( F, norm( reshape( (ARay_after_calibration.*repmat( AClosestDepth, 3,1)) - ...
                           (R2'*(BRay_after_calibration.*repmat( BClosestDepth, 3,1)) + repmat(T_simple_before_R, 1, length(ADepthMatch))) , 1,[]), 1), opt);
    AClosestDepth = double(AClosestDepth);
    BClosestDepth = double(BClosestDepth);
    
    
    % check if the ray really match well =============================
    figure(50); title('Closest point Match Point, one opt with BScale');
    AClosestMatchPosition = ARay_after_calibration.*repmat( AClosestDepth, 3,1);
    BClosestMatchPosition = R2'*(BRay_after_calibration.*repmat( BClosestDepth, 3,1)) + repmat(T_simple_before_R, 1, length(ADepthMatch));
    scatter3(AClosestMatchPosition(1,:)', AClosestMatchPosition(3,:)', AClosestMatchPosition(2,:)', 0.5*ones(1,size( AClosestMatchPosition,2)));
    hold on;
    scatter3(BClosestMatchPosition(1,:)', BClosestMatchPosition(3,:)', BClosestMatchPosition(2,:)', 0.5*ones(1,size( BClosestMatchPosition,2)));
    line( [ AClosestMatchPosition(1,:); BClosestMatchPosition(1,:)], ...
          [ AClosestMatchPosition(3,:); BClosestMatchPosition(3,:)], ...
          [ AClosestMatchPosition(2,:); BClosestMatchPosition(2,:)]);
    % ================================================================
 end
    % check if the ray really match well =============================
    figure(52); title('Match Point, one opt with BScale');
    ADepth_modifiedPosition = ARay_after_calibration.*repmat( ADepth_modified, 3,1);
    BDepth_modifiedPosition = R2'*(BRay_after_calibration.*repmat( BDepth_modified, 3,1)) + repmat(T_simple_before_R, 1, length(ADepthMatch));
    scatter3(ADepth_modifiedPosition(1,:)', ADepth_modifiedPosition(3,:)', ADepth_modifiedPosition(2,:)', 0.5*ones(1,size( ADepth_modifiedPosition,2)));
    hold on;
    scatter3(BDepth_modifiedPosition(1,:)', BDepth_modifiedPosition(3,:)', BDepth_modifiedPosition(2,:)', 0.5*ones(1,size( BDepth_modifiedPosition,2)));
    line( [ ADepth_modifiedPosition(1,:); BDepth_modifiedPosition(1,:)], ...
          [ ADepth_modifiedPosition(3,:); BDepth_modifiedPosition(3,:)], ...
          [ ADepth_modifiedPosition(2,:); BDepth_modifiedPosition(2,:)]);    
    % ================================================================

    BTranslation = ATranslation + ARotation*T_simple_before_R;
% =========not used anymore since scaling is included in match point jointly opt
%     % solve the scale for depth A and depth B 
%     ADepthScale = sdpvar(1); 
%     BDepthScale = sdpvar(1); 
%     sol =solvesdp( [], norm( ADepth_modified - ADepthScale*ADepthMatch, 1)+...
%                     norm( BDepth_modified - BDepthScale*BDepthMatch, 1),opt);
%     ADepthScale = double(ADepthScale);
%     BDepthScale = double(BDepthScale);
%     %ADepth_normalized = ADepth*ADepthScale;
%     %BDepth_normalized = BDepthMatch*BDepthScale;
% =========================================================================
    
    if RenderWrlDirectlyFlag % only a reality check to see if the rotation and translation is reasonable
         % rendering the joint wrl
         temp = Ray(:,:,1:2)./repmat(Ray(:,:,3),[1 1 2]);
         PositionTex = permute(temp./repmat(cat(3,Default.a_default,Default.b_default),[Default.VertYNuDepth Default.HoriXNuDepth 1])+repmat(cat(3,Default.Ox_default,Default.Oy_default),[Default.VertYNuDepth Default.HoriXNuDepth 1]),[3 1 2]);
         PositionTex = permute(PositionTex,[2 3 1]);
         %  if SeperateOptFlag
         AWrlPosition = ARotation*APositionAll(:,:)*ADepthScale+...
                        repmat( ATranslation, 1, size(APositionAll(:,:), 2));
         BWrlPosition = BRotation*BPositionAll(:,:)*BDepthScale+...
                        repmat( BTranslation, 1, size(BPositionAll(:,:),2));
         AWrlPosition = reshape(AWrlPosition,3,55,[]);
         AWrlPosition(3,:) = - AWrlPosition(3,:);
         AWrlPosition = permute(AWrlPosition,[2 3 1]);
         BWrlPosition = reshape(BWrlPosition,3,55,[]);
         BWrlPosition(3,:) = - BWrlPosition(3,:);
         BWrlPosition = permute(BWrlPosition,[ 2 3 1]);
         WrlFacestHroiReduce(AWrlPosition,PositionTex,ASup,left,[ Wrlname '_' left '-' right '_OldMulti'], [ResultFolder '/'], 0,0);
         WrlFacestHroiReduce(BWrlPosition,PositionTex,BSup,right,[ Wrlname '_' left '-' right '_OldMulti'], [ResultFolder '/'], 0,1);
         system(['cp ' Fdir '/jpg/' left '.jpg ' ResultFolder '/' left '.jpg']);
         system(['cp ' Fdir '/jpg/' right '.jpg ' ResultFolder '/' right '.jpg']);
         disp('Finish model simple rotation and translation, first reality check');
%         pause;
    end
    
% =========estimating ground level using ground in imgA and imgB on image A  coordinate
 % cleaning Groundmask for A
RangePercent = 100;
AY_median = median( APositionAll(2,AmaskG));
ADistant2Ay_median = (APositionAll(2,AmaskG) - AY_median);
ANumber_YMedian = round( sum(APositionAll(2,AmaskG)<AY_median)*RangePercent/100);
[Avalue AIndexSort] = sort(ADistant2Ay_median);
AYmedia_mark = AIndexSort(1:ANumber_YMedian);
AGround_mark = zeros(Dy, Dx);
temp = zeros(sum(AmaskG(:)) ,1);
temp(AYmedia_mark) = 1;
AGround_mark(AmaskG) = temp;
AGround_mark = logical(AGround_mark);
   % finishing cleaning Groundmask for A
   
 % cleaning Groundmask for B
BWrlPosition = R2'*BPositionAll(:,:)*BDepthScale+repmat( T_simple_before_R, 1, size(BPositionAll(:,:),2)); 
BWrlPosition = reshape(BWrlPosition,3,55,[]);
BY_median = median( BWrlPosition(2,BmaskG));
BDistant2By_median = (BWrlPosition(2,BmaskG) - BY_median);
BNumber_YMedian = round( sum(BWrlPosition(2,BmaskG)<BY_median)*RangePercent/100);
[Bvalue BIndexSort] = sort(BDistant2By_median);
BYmedia_mark = BIndexSort(1:BNumber_YMedian);
BGround_mark = zeros(Dy, Dx);
temp = zeros(sum(BmaskG(:)) ,1);
temp(BYmedia_mark) = 1;
BGround_mark(BmaskG) = temp;
BGround_mark = logical( BGround_mark);
 % finishing cleaning Groundmask for B
 
% find the jointly median of the ground of image AB in Y direction
if ~appendOpt
   GroundY0inACoor = median([ APositionAll(2,AGround_mark) BWrlPosition(2,BGround_mark)]);
end
% =====================================

    % given the AClosestDepth and BClosestDepth do the seperate inference again output the ADepth_normalized and BDepth_normalized
    % imgA
    AappendOpt = appendOpt;
    ASupMatched = ASup(Aindex)';
    mask = ASupMatched == 0;
    ASupMatched(mask)=[];
    ARayMatched = ARay_after_calibration';
    ARayMatched(mask,:) = [];
    ADepth_modified(:,mask) = [];
%     AClosestDepth(:,mask) = [];
    Default.OutPutFolder = OutPutFolder;
    Default.ScratchFolder = ScratchFolder;
    Default.filename{1} = left;
    Default.Wrlname{1} = Wrlname;
    Default.Flag.AfterInferenceStorage = 0;
%    if ~appendOpt
    if RenderFlag
       Default.RenderFlag = 1;
    else
       Default.RenderFlag = 0;
    end
     [APosition_normalized ADepth_normalized] = ReportPlaneParaMRF_Conditioned_trianglate2( Default, ARotation, ATranslation, AappendOpt, ...
                           [ ARayMatched; Aconstrain.RayMatched],...
                           [ ADepth_modified Aconstrain.Depth_modified]',...
                           [ ASupMatched; Aconstrain.SupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, AMedSup, AdepthMap*ADepthScale, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, 0, ADirectionFromReference, GroundY0inACoor);
%    else
%        ADepth_normalized = AdepthMap;
%    end
%    [APosition_normalized ADepth_normalized] = ReportPlaneParaMRF_Conditioned_trianglate( Default, ARotation, ATranslation, AappendOpt, ...
%                          [ ], [ ], [ ], ...
%                          [ ], [ ], [ ], [],...
%                          ASup, ASupOri, AMedSup, AdepthMap*ADepthScale, ones(size(AdepthMap)), ARayOri, ARayAll, ...
%                          ASupNeighborTable, [], AmaskSky, AmaskG,...
%                          'cvx_allL1Norm',1,...
%                          [], [], AMultiScaleSupTable, [], [], [], false, 0, [], GroundY0inACoor);
    % imgB
%    BRotation = R2'*ARotation';
%    BTranslation = ATranslation + ARotation*T_simple_before_R;
    BappendOpt = 1;
    BSupMatched = BSup(Bindex)';
    mask = BSupMatched == 0;
    BSupMatched(mask)=[];
    BRayMatched = BRay_after_calibration';
    BRayMatched(mask,:) = [];
    BDepth_modified(:,mask) = [];
%     BClosestDepth(:,mask) = [];
    Default.OutPutFolder = OutPutFolder;
    Default.ScratchFolder = ScratchFolder;
    Default.filename{1} = right;
    Default.Wrlname{1} = Wrlname;
    Default.Flag.AfterInferenceStorage = 0;
    
    Default.RenderFlag = 0;
     [BPosition_normalized BDepth_normalized] = ...
                           ReportPlaneParaMRF_Conditioned_trianglate2( Default, BRotation, BTranslation, BappendOpt, ...
                           [ BRayMatched ], [ BDepth_modified' ], [ BSupMatched], ...
                           [ ], [ ], [ ], [],...
                           BSup, BSupOri, BMedSup, BdepthMap*BDepthScale, zeros(size(BdepthMap)), BRayOri, BRayAll, ...
                           BSupNeighborTable, [], BmaskSky, BmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], BMultiScaleSupTable, [], [], [], false, 0, BDirectionFromReference, GroundY0inACoor);
%   [BPosition_normalized BDepth_normalized] = ...
%                          ReportPlaneParaMRF_Conditioned_trianglate( Default, BRotation, BTranslation, BappendOpt, ...
%                          [  ], [  ], [ ], ...
%                          [  ], [  ], [ ], [],...
%                          BSup, BSupOri, BMedSup, BdepthMap*BDepthScale, ones(size(BdepthMap)), BRayOri, BRayAll, ...
%                          BSupNeighborTable, [], BmaskSky, BmaskG,...
%                          'cvx_allL1Norm',1,...
%                          [], [], BMultiScaleSupTable, [], [], [], false, 0, BDirectionFromReference, GroundY0inACoor);                      
%  
    AMatchedPosition_after_normalized = (ARay_after_calibration').*repmat( ADepth_normalized(Aindex)', 1, 3);
    BMatchedPosition_after_normalized = (BRay_after_calibration').*repmat( BDepth_normalized(Bindex)', 1, 3);

% ============= save constrain Rotation Translation and combining history ===============
%ImgA
AHistory{end+1}= right;
if isempty(intersect(AHistory, right))
  Aconstrain.RayMatched = [ Aconstrain.RayMatched; ARayMatched ];
  Aconstrain.Depth_modified = [ Aconstrain.Depth_modified ADepth_modified ];
  Aconstrain.SupMatched = [ Aconstrain.SupMatched; ASupMatched ];
end  
AdepthMap = ADepth_normalized;
save([ ScratchFolder '/' left '_NonMono.mat' ], 'AdepthMap', 'ARotation', 'ATranslation', 'AHistory', 'Aconstrain',...
             'GroundY0inACoor',...
            'ASup', 'ASupOri', 'AMedSup', 'ARayOri','ARayAll','ASupNeighborTable','AmaskSky','AmaskG','AMultiScaleSupTable');

%ImgB
ARotation = BRotation;
ATranslation = BTranslation;
AHistory{1}= 'left';
AdepthMap = BdepthMap;
ARayAll = BRayAll;
ARayOri = BRayOri;
ASup = BSup;
ASupOri = BSupOri;
AMedSup = BMedSup;
ASupNeighborTable = BSupNeighborTable;
AmaskSky = BmaskSky;
AmaskG = BmaskG;
AMultiScaleSupTable = BMultiScaleSupTable;
Aconstrain.RayMatched = [ BRayMatched ];
Aconstrain.Depth_modified = [ BDepth_modified ];
Aconstrain.SupMatched = [ BSupMatched ];
AdepthMap = BDepth_normalized;
save([ ScratchFolder '/' right '_NonMono.mat' ], 'AdepthMap', 'ARotation', 'ATranslation', 'AHistory', 'Aconstrain',...
             'GroundY0inACoor',...
            'ASup', 'ASupOri', 'AMedSup', 'ARayOri','ARayAll','ASupNeighborTable','AmaskSky','AmaskG','AMultiScaleSupTable');


% =======================================================================================
 % generate the 3D model jointly
% Aaxis = [ zeros(1,3); [0 0 10]; zeros(1,3); [0 10 0]];
% BCenter = T_simple_before_R;
% Baxis = [ BCenter'; (R2'*[0 0 10]'+BCenter)'; BCenter'; (R2'*[0 10 0]'+BCenter)'];
%  DisplayPairPointsCloud(APositionAll(:,:)', (R2'*BPositionAll(:,:)+repmat( T_simple_before_R, 1, size(BPositionAll(:,:),2)))',...
%                        Aaxis, Baxis, 300  );
%  DisplayPairPointsCloud(APositionAll(:,:)'*ADepthScale, (R2'*BPositionAll(:,:)*BDepthScale+repmat( T_simple_before_R, 1, size(BPositionAll(:,:),2)))',...
%                        Aaxis, Baxis, 400  );                  
                   

 
%Psimple = Pnew*Hsimple;
%Xsimple = inv(Hsimple)*Xnew;
%DisplayMetricReconstruction( Psimple,Xsimple);
%DisplayMetricReconstruction( Psimple,Xsimple);
return;
% ========================================================
