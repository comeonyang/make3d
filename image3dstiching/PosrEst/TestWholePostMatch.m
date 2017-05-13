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
function [Pair fail]=TestWholePostMatch(defaultPara, ImgInfo1, ImgInfo2, PriorPose)

% This function load the Initial SurfMatches and the OccluSurfMatches
% then, do RANSAC algorithm to prune the bad matches out
% Calucluate New Camera Pose (R and T)

% Default Parameters -----------------------------
fail = 0;
% ------------------------------------------------

% Initialize Variables ---------------------------
Img1 = ImgInfo1.ExifInfo.IDName;
Img2 = ImgInfo2.ExifInfo.IDName;
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
% ------------------------------------------------

% Combine Both Initial SurfMatches and the OccluSurfMatches -------------------------
% Also Assign Confidence of each Matches 
% (Used as RANSAC Sampling Distribution)
% 0) load cleaned Fisrt set of matches
%defaultPara.MaxUniqueRatio
load([defaultPara.Fdir '/data/' Img1 '_' Img2 '_PoseMatch.mat']);
defaultPara.MaxUniqueRatio = 100;
% 1) load Initial SurfMatches
[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, ...
                    [ defaultPara.Type 'Dense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)], 1, 1);
if isempty(f1)
	[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, ...
                        ['Dense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)], 1, 1);
end

% 2) load OccluSurfMatches
[f1, f2, OccluedMatches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, [ defaultPara.Type 'OccluDense'], 1, 1, 3);% need UniqueRatio//Min jobs
UniqueRatio = OccluedMatches(3,:);
%UniqueRatio = min(UniqueRatio,defaultPara.MaxUniqueRatio);
Ptr = UniqueRatio > defaultPara.MaxUniqueRatio;
UniqueRatio(Ptr) = Inf;%defaultPara.MaxUniqueRatio;

OccluedMatches = OccluedMatches(1:2,:);
if false
	disp('use CleanMatches');
	matches = Matches;
	% removing initial matches that inconsistent with OccluedMatches
	[c i] = intersect( matches(1,:), OccluedMatches(1,:));
	matches(:,i) = [];
	[c i] = intersect( matches(2,:), OccluedMatches(2,:));
	matches(:,i) = [];
	matches = [matches OccluedMatches];

else
	disp('Number of Occlumatches used')
	size( OccluedMatches,2)

	% removing initial matches that inconsistent with OccluedMatches
	[c i] = intersect( matches(1,:), OccluedMatches(1,:));
	matches(:,i) = [];
	[c i] = intersect( matches(2,:), OccluedMatches(2,:));
	matches(:,i) = [];

	NumInitialMatches = size(matches, 2);
	matches = [matches OccluedMatches];
	if isempty(matches)
	   disp('Zeros Surf matches');
	   fail = 1;
	   return;
	end    

	% 3) Construct Prior Dist for All matches
	NMatches = size(matches, 2);
	PriorDist = ones(1, NMatches);
	PriorDist( (NumInitialMatches+1):end) = UniqueRatio;
	% 4) RANSAC 
	disp('Number of matches used')
	
	if defaultPara.Flag.StorageDataBeforeRansac
		disp('Storage Data Before Ransac');
		save([defaultPara.Fdir '/data/DataBeforeRansac.mat']);
% 		return;
	end

	% Ensemble method to determine confidence of inliers
	fittingfn = @fundmatrix;
	distfnEnsmble = @fundistEnsmble;
	degenfn   = @isdegenerate;
	x = [[f1(:, matches(1,:)); ones(1, NMatches)]; [f2(:, matches(2,:)); ones(1, NMatches)]];
	[ SampsonDist ] = EnsembleRansac(defaultPara, x, fittingfn, distfnEnsmble, degenfn, 8, PriorDist', min(NMatches*10, defaultPara.MAXEnsembleSamples), 0);
	kurtosisValue =kurtosis(SampsonDist');

	[F0, inliers, NewDist, fail, ind]=GeneralRansac(defaultPara, f1, f2, matches, [], [], kurtosisValue', 8);
	figure(100); plotmatches(I1,I2,f1, f2,matches(:,inliers), 'Stacking', 'h', 'Interactive', 3);
	matches = matches(:,inliers); % accept the pruning result
	if isempty(matches)
	   disp('Zeros After Ransac matches');
	   fail = 2;
	   return;
	end
end
% ------------------------------------------------------------------------------------


% Apply Bundle Adjustment Refinment Algorithm to Prune the Matches further -----------
% And Estimated the Pose
% 1) initialize the 3D position of the matches given Prior Pose and Depths
x_calib = [ inv(defaultPara.InrinsicK1)*[ f1(:,matches(1,:)); ones(1,size(matches,2))];...
	      inv(defaultPara.InrinsicK2)*[ f2(:,matches(2,:)); ones(1,size(matches,2))]];

% Estimate F using NonLine LS on every inlier
MatchDensityWeights1 = CalMatchDensityWeights(f1(:,matches(1,:)), max(size(I1))/defaultPara.radius2imageSizeRatio);
MatchDensityWeights2 = CalMatchDensityWeights(f2(:,matches(2,:)), max(size(I2))/defaultPara.radius2imageSizeRatio);
MatchDensityWeights =mean([MatchDensityWeights1; MatchDensityWeights2], 1);
F = getFnpt( F0, f1(:, matches(1,:))',  f2(:, matches(2,:))', MatchDensityWeights);
E = defaultPara.InrinsicK2'*F*defaultPara.InrinsicK1; % Camera essential Matrix
if ~isempty(PriorPose)
	[ R0, T0, lamda1, lamda2, inlier, Error] = EstPose( defaultPara, E, x_calib, [], PriorPose.R(1:3,:));
else
	[ R0, T0, lamda1, lamda2, inlier, Error] = EstPose( defaultPara, E, x_calib, [], []);
end
T0 = [T0; - R0'*T0];
R0 = [R0; R0'];
matches = matches(:,inlier); % delet matches give negative depths
x_calib = x_calib(:,inlier);
lamda1 = lamda1(inlier);
lamda2 = lamda2(inlier);

X_obj_1 = x_calib(1:3,:).*repmat(lamda1, 3, 1);
X_obj_2 = R0(4:6,:)*(x_calib(4:6,:).*repmat(lamda2, 3, 1)) + repmat(T0(4:6), 1, size(matches,2));
X_obj = (X_obj_1+X_obj_2)/2;

% 2) 
[R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, R0(1:3,:), T0(1:3), [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj, [ ImgInfo1 ImgInfo2], 1);
if false % Min Modified for not pruning using BundleAdjustment
	if all(isnan( dist1_BA)) || isempty(R) || any(isnan(R(:)))
		disp('BA failed');
		fail = 3;
		return;
	end
	while length(X_im_BA) >= defaultPara.MinimumNumMatches
		disp('BundleAdjClean')
	    outlier_thre1 = prctile(dist1_BA,90);
	    outlier_thre2 = prctile(dist2_BA,90);
	    if outlier_thre1 >= defaultPara.ReProjErrorThre || outlier_thre2 >= defaultPara.ReProjErrorThre
		Outlier = dist1_BA > outlier_thre1 | dist2_BA > outlier_thre2;
		matches(:,Outlier) = [];
		if isempty(matches)
		    disp('Zeros After BA matches');
		    fail = 4;
		    return;
		end
		lamda1(Outlier) = [];
		lamda2(Outlier) = [];
		X_obj_BA(:,Outlier) = [];
		x_calib(:,Outlier) = [];
		[R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, R, T, [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj_BA, [ ImgInfo1 ImgInfo2], 1);
		if all(isnan( dist1_BA)) || isempty(R) || any(isnan(R(:)))
			disp('BA failed');
		    fail = 5;
			return;
		end
	    else
		break;
	    end
	end    
end
% ------------------------------------------------------------------------------------


% Triangulation ----------------------------------------------------------------------
% modified the x_calib
tempf1 = X_im_BA(1:2,:);
tempf2 = X_im_BA(3:4,:);
x_calib = [ inv(defaultPara.InrinsicK1)*[ tempf1; ones(1,length(tempf1))];...
             inv(defaultPara.InrinsicK2)*[ tempf2; ones(1,length(tempf2))]];
[ lamda1 lamda2 Error] = triangulation( defaultPara, R, T, x_calib);
% Clean outlier triangulated depth
%LamdaOutlier = lamda1 > 1000 | lamda1 <1;
%matches(:,LamdaOutlier) = [];
%if isempty(matches)
%    disp('Zeros After Tri matches');
%    fail = 6;
%    return;
%end

%lamda1(LamdaOutlier) = [];
%lamda2(LamdaOutlier) = [];
%x_calib(:,LamdaOutlier) = [];

% Pair Image Depth Scale
[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo1.Model.Depth.FitDepth), f1(:,matches(1,:)), ImgInfo1.Model.Depth.FitDepth);
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo2.Model.Depth.FitDepth), f2(:,matches(2,:)), ImgInfo2.Model.Depth.FitDepth);
Depth1ProjDepthRatio = sqrt(sum(x_calib(1:3,:).^2, 1));
Depth2ProjDepthRatio = sqrt(sum(x_calib(4:6,:).^2, 1));
DProj1 = D1./Depth1ProjDepthRatio;
DProj2 = D2./Depth2ProjDepthRatio;
[DepthScale1] = UniformDepthScale( defaultPara, DProj1, lamda1, ones(1,length(lamda1)));
[DepthScale2] = UniformDepthScale( defaultPara, DProj2, lamda2, ones(1,length(lamda2)) );
%if DepthScale1 > 20 | DepthScale1 <0.05 | DepthScale2 > 20 | DepthScale2 <0.05 %//Min used to use 10 and 0.2
%   disp('Unrealistic in Rescaleing the depth, Check matchings');
%   fail = 6;
%   return;
%end 

Pair.lamda = [lamda1; lamda2];%//Min add for debug
Pair.DepthScale = [DepthScale1; DepthScale2];
Pair.R = R;
Pair.T = T;
Pair.Xim = [f1(:,matches(1,:)); f2(:,matches(2,:))];
% -----------------------------------------------------------------------------------
