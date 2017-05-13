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
function [Pair ImgInfo matches fail]=PoseMatchEst(defaultPara, ImgInfo)

% This function estimate the relative Pose of the camera using first camera coordinate
% as world coordinate

% Input:
% 	default - camera intrinsic, etc
%	ImgInfo - Exif, Model info, GPS, IMU info
%
% Return:	
%	R - rotation - (R*Posi2+ T to A's coordinate)
%	T - translation

% step outline
%	1) extract Measuesd Position and orientation from GPS or IMU info
%	2) Using Measures R and T and Mono-Depth to define mach search space constrain
%	3) Do match search with all combinations satisfying Constrain from 2) using ralative threshould
%	4) Ransac
%	5) Bundle Adjustment
%	6) up to scale translation reconstruction
%	7) matches 3D triangulation
%	8) Modified ImgInfo.Model.Depth up to accurate scale

% initialize variables
fail = 0;
Pair.R = [];
Pair.t = [];
Pair.Xim = [];
Pair.DepthScale = [];

Img1 = ImgInfo(1).ExifInfo.IDName;
Img2 = ImgInfo(2).ExifInfo.IDName;
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
[f1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % Dense features 
[f2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % Dense features 
[D1 IND] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1, ImgInfo(1).Model.Depth.FitDepth);
[D2 IND] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2, ImgInfo(1).Model.Depth.FitDepth);

% 1) extract Measuesd Position and orientation from GPS or IMU info
% Depends on what data we have, MeasR or MeasT, or both might be empty
[MeasR MeasT] = InitPoseMeas(defaultPara, ImgInfo(1), ImgInfo(2));

if ~isempty(MeasR)
	% 2) Using Measures R and T and Mono-Depth to define match search space constrain
	% read in all surf features
	[ Rc1, Rc2, ConS1, ConS2, ConSRough1, ConSRough2] = CalMatchSearchRegin(defaultPara, MeasR, MeasT, I1, I2, f1, f2, D1, D2, 1, defaultPara.Flag.FlagDisp);

	% write the match search space constrain in to files for surfMatchRConS.sh script to read
	Vector2Ipoint([Rc1; ConS1],[defaultPara.Fdir '/surf/'],['RConS_' Img1]);
	Vector2Ipoint([Rc2; ConS2],[defaultPara.Fdir '/surf/'],['RConS_' Img2]);
	Vector2Ipoint([ConSRough1],[defaultPara.Fdir '/surf/'],['RConSRough_' Img1]);
	Vector2Ipoint([ConSRough2],[defaultPara.Fdir '/surf/'],['RConSRough_' Img2]);

	% 3) Do match search with all combinations satisfying Constrain from 2) using ralative threshould
	cd match
	[status, result] = system(['ls ' defaultPara.Fdir '/surf_matches/' Img1 '-' Img2 '.matchRConSDense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)]);
	[statusReverse, resultReverse] ...
			 = system(['ls ' defaultPara.Fdir '/surf_matches/' Img2 '-' Img1 '.matchRConSDense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)]);
	if status && statusReverse
		SurfMatchTime = tic;
		system(['./surfMatchRConS.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' Dense ' num2str(defaultPara.AbsThre) ' ' num2str(defaultPara.RatioThre)]); 
		disp(['		' num2str( toc( SurfMatchTime)) ' seconds.']);
	end    
	cd ..
else
	cd match
	[status, result] = system(['ls ' defaultPara.Fdir '/surf_matches/' Img1 '-' Img2 '.matchDense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)]);
	[statusReverse, resultReverse] ...
			 = system(['ls ' defaultPara.Fdir '/surf_matches/' Img2 '-' Img1 '.matchDense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)]);
	if status && statusReverse
		SurfMatchTime = tic;
		system(['./surfMatch.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' Dense ' num2str(defaultPara.AbsThre) ' ' num2str(defaultPara.RatioThre)]); 
		disp(['		' num2str( toc( SurfMatchTime)) ' seconds.']);
	end    
	cd ..
end
% 4. Ransac
if ~isempty(MeasR)
	[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, ...
				[ defaultPara.Type 'Dense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)], 1, 1);
else
	[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, ...
				[ 'Dense_' num2str(defaultPara.AbsThre) '_' num2str(defaultPara.RatioThre)], 1, 1);
end
if isempty(matches)
   disp('Zeros Surf matches');
   fail = 1;
   return;
end    
[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1(:,matches(1,:)), ImgInfo(1).Model.Depth.FitDepth);
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2(:,matches(2,:)), ImgInfo(1).Model.Depth.FitDepth);

% Ensemble method to determine confidence of inliers
fittingfn = @fundmatrix;
distfnEnsmble = @fundistEnsmble;
degenfn   = @isdegenerate;
nmatches = size(matches, 2);
x = [[f1(:, matches(1,:)); ones(1, nmatches)]; [f2(:, matches(2,:)); ones(1, nmatches)]];
[ SampsonDist ] = EnsembleRansac(defaultPara, x, fittingfn, distfnEnsmble, degenfn, 8, ones(1,nmatches)', min(nmatches*10, defaultPara.MAXEnsembleSamples), 0);
kurtosisValue =kurtosis(SampsonDist');

% Ransac
[F0, inliers, NewDist, fail, ind]=GeneralRansac(defaultPara, f1, f2, matches, D1, D2, kurtosisValue', 8);
if defaultPara.Flag.FlagDisp
    figure; plotmatches(I1,I2,f1, f2,matches(:,inliers), 'Stacking', 'v', 'Interactive', defaultPara.Flag.FlagDisp);
end   
% *** Stop maunally to pick out the bad matches*** -----------------
matches = matches(:,inliers);
if isempty(matches)
   disp('Zeros Matches After Ransac');
   fail = 2;
   return;
end
% ------------------------------------------------------------------

x_calib = [ inv(defaultPara.InrinsicK1)*[ f1(:,matches(1,:)); ones(1,length(matches))];...
	      inv(defaultPara.InrinsicK2)*[ f2(:,matches(2,:)); ones(1,length(matches))]]; 

% Estimate F using NonLine LS on every inlier
MatchDensityWeights1 = CalMatchDensityWeights(f1(:,matches(1,:)), max(size(I1))/defaultPara.radius2imageSizeRatio);
MatchDensityWeights2 = CalMatchDensityWeights(f2(:,matches(2,:)), max(size(I2))/defaultPara.radius2imageSizeRatio);
MatchDensityWeights =mean([MatchDensityWeights1; MatchDensityWeights2], 1);
F = getFnpt( F0, f1(:, matches(1,:))',  f2(:, matches(2,:))', MatchDensityWeights);
E = defaultPara.InrinsicK2'*F*defaultPara.InrinsicK1; % Camera essential Matrix
if ~isempty(MeasR)
	[ R0, T0, lamda1, lamda2, inlier, Error] = EstPose( defaultPara, E, x_calib, [], MeasR(1:3,:));
else
	[ R0, T0, lamda1, lamda2, inlier, Error] = EstPose( defaultPara, E, x_calib, [], []);
end
T0 = [T0; - R0'*T0];
R0 = [R0; R0'];
matches = matches(:,inlier); % delet matches give negative depths
x_calib = x_calib(:,inlier);
lamda1 = lamda1(inlier);
lamda2 = lamda2(inlier);

% Estimated X_obj by triangulation
X_obj_1 = x_calib(1:3,:).*repmat(lamda1, 3, 1);
X_obj_2 = R0(4:6,:)*(x_calib(4:6,:).*repmat(lamda2, 3, 1)) + repmat(T0(4:6), 1, size(matches,2));
X_obj = (X_obj_1+X_obj_2)/2;

% 5. Bundle Adjustment
[R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, R0(1:3,:), T0(1:3), [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj, ImgInfo, 1);
if all(isnan( dist1_BA)) || isempty(R) || any(isnan(R(:)))
	disp('BA failed');
   	fail = 3;
        return;
end
while length(X_im_BA) >= defaultPara.MinimumNumMatches
	outlier_thre1 = prctile(dist1_BA,90);
	outlier_thre2 = prctile(dist2_BA,90);
	Outlier = logical(zeros( size( dist1_BA)));
	if max(dist1_BA) >= defaultPara.ReProjErrorThre 
% 		Outlier = Outlier | dist1_BA > max( outlier_thre1, defaultPara.ReProjErrorThre);
		Outlier = Outlier | dist1_BA > outlier_thre1;        
	end
	if max(dist2_BA) >= defaultPara.ReProjErrorThre
% 		Outlier = Outlier | dist2_BA > max( outlier_thre2, defaultPara.ReProjErrorThre);
        Outlier = Outlier | dist2_BA > outlier_thre2;
	end
	matches(:,Outlier) = [];
	if isempty(matches)
	    disp('Zeros Matches After BA Pruning');
	    fail = 4;
	    return;
	end
	if all( Outlier == 0)
		% Non Outlier detected for BA
		break;
	end
	lamda1(Outlier) = [];
	lamda2(Outlier) = [];
	X_obj_BA(:,Outlier) = [];
	x_calib(:,Outlier) = [];
	[R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, R, T, [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj_BA, ImgInfo, 1);
	if all(isnan( dist1_BA)) || isempty(R) || any(isnan(R(:)))
		disp('BA failed');
	    fail = 5;
		return;
	end
end    
if defaultPara.Flag.FlagDisp
    figure;  plotmatches(I1,I2,f1, f2, matches, 'Stacking', 'v', 'Interactive', defaultPara.Flag.FlagDisp);
end
    
% 6. find T up to scale

% 7. Triangulation 
% modified the x_calib So that perfact triangulation but the image is distorted a little bit
tempf1 = X_im_BA(1:2,:);
tempf2 = X_im_BA(3:4,:);
x_calib = [ inv(defaultPara.InrinsicK1)*[ tempf1; ones(1,length(tempf1))];...
             inv(defaultPara.InrinsicK2)*[ tempf2; ones(1,length(tempf2))]];
% ------------------
[ lamda1 lamda2 Error] = triangulation( defaultPara, R, T, x_calib);

% 8. modify ImgInfo.Model.Depth .... (not sure do it or not??????)
[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1(:,matches(1,:)), ImgInfo(1).Model.Depth.FitDepth);
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2(:,matches(2,:)), ImgInfo(2).Model.Depth.FitDepth);
Depth1ProjDepthRatio = sqrt(sum(x_calib(1:3,:).^2, 1));
Depth2ProjDepthRatio = sqrt(sum(x_calib(4:6,:).^2, 1));
DProj1 = D1./Depth1ProjDepthRatio;
DProj2 = D2./Depth2ProjDepthRatio;
[DepthScale1] = UniformDepthScale( defaultPara, DProj1, lamda1, ones(1,length(lamda1)));
[DepthScale2] = UniformDepthScale( defaultPara, DProj2, lamda2, ones(1,length(lamda2)) );
%if DepthScale1 > 20 | DepthScale1 <0.05 | DepthScale2 > 20 | DepthScale2 <0.05 %//Min used to use 10 and 0.2
%   disp('Unrealistic in Rescaleing the depth, Check matchings');
%   fail = -1;
%end 

Pair.lamda = [lamda1; lamda2];
Pair.DepthScale = [DepthScale1; DepthScale2];
Pair.R = R;
Pair.T = T;
Pair.Xim = [f1(:,matches(1,:)); f2(:,matches(2,:))];

% check is triangulation reasonable
if defaultPara.Flag.FlagDisp
figure(50); clf; title('Closest point Match Point'); hold on;
ClosestMatchPosition2 = x_calib(4:6,:).*repmat( lamda2, 3,1);
ClosestMatchPosition1 = R*(x_calib(1:3,:).*repmat( lamda1, 3,1)) + repmat(T, 1, length(lamda1));
MonoStichPosition2 = x_calib(4:6,:).*repmat( DProj2.*DepthScale2, 3,1);
MonoStichPosition1 = R*(x_calib(1:3,:).*repmat( DProj1.*DepthScale1, 3,1)) + repmat(T, 1, length(DProj1));
% =====================
[VDepth HDepth] = size(ImgInfo(2).Model.Depth.FitDepth);
[VImg HImg] = size(I1);
VIndexDepthRes = repmat((1:VDepth)', [1 HDepth]);
HIndexDepthRes = repmat((1:HDepth), [VDepth 1]);
VIndexImgRes = ( VIndexDepthRes -0.5)/VDepth*VImg;
HIndexImgRes = ( HIndexDepthRes -0.5)/HDepth*HImg;
ImgPositionPix = cat(3, HIndexImgRes, VIndexImgRes);
All_x_calib = inv(defaultPara.InrinsicK1)*[ reshape( permute(ImgPositionPix, [ 3 1 2]), 2, []); ones(1, VDepth*HDepth)];%
All_Ray = All_x_calib./repmat( sqrt(sum(All_x_calib.^2, 1)), 3, 1);
All_Ray = repmat( All_Ray, 2, 1);
% ====================
ReScaledPosi2 = All_Ray(4:6,:).*repmat( ImgInfo(2).Model.Depth.FitDepth(:)'*DepthScale2, 3,1);
ReScaledPosi1 = R*(All_Ray(1:3,:).*repmat( ImgInfo(1).Model.Depth.FitDepth(:)'*DepthScale1, 3,1)) + repmat(T, 1, length(All_Ray));
ReScaledPosi2(:,IND2) = [];
ReScaledPosi1(:,IND1) = [];
scatter3(ReScaledPosi2(1,:)', ReScaledPosi2(3,:)', ReScaledPosi2(2,:)', 0.5*ones(1,size( ReScaledPosi2,2)));
scatter3(ReScaledPosi1(1,:)', ReScaledPosi1(3,:)', ReScaledPosi1(2,:)', 1*ones(1,size( ReScaledPosi1,2)));
scatter3(ClosestMatchPosition2(1,:)', ClosestMatchPosition2(3,:)', ClosestMatchPosition2(2,:)', 40, 'g');
scatter3(ClosestMatchPosition1(1,:)', ClosestMatchPosition1(3,:)', ClosestMatchPosition1(2,:)', 40, 'b');
line( [ ClosestMatchPosition2(1,:); ClosestMatchPosition1(1,:)], ...
      [ ClosestMatchPosition2(3,:); ClosestMatchPosition1(3,:)], ...
      [ ClosestMatchPosition2(2,:); ClosestMatchPosition1(2,:)]);
% line( [ MonoStichPosition2(1,:); MonoStichPosition1(1,:)], ...
%       [ MonoStichPosition2(3,:); MonoStichPosition1(3,:)], ...
%       [ MonoStichPosition2(2,:); MonoStichPosition1(2,:)]);
if ~isempty(ImgInfo(1).Model.Constrain.RayMatche)
  ClosestMatchPosition1Hist = R*(ImgInfo(1).Model.Constrain.RayMatche'.*repmat(ImgInfo(1).Model.Constrain.Depth_modified , 3, 1)) + repmat(T, 1, length(ImgInfo(1).Model.Constrain.RayMatche));
  scatter3(ClosestMatchPosition1Hist(1,:)', ClosestMatchPosition1Hist(3,:)', ClosestMatchPosition1Hist(2,:)', 40, 'y');
end
if ~isempty(ImgInfo(2).Model.Constrain.RayMatche)
  ClosestMatchPosition2Hist = ImgInfo(2).Model.Constrain.RayMatche'.*repmat(ImgInfo(2).Model.Constrain.Depth_modified , 3, 1);
  scatter3(ClosestMatchPosition2Hist(1,:)', ClosestMatchPosition2Hist(3,:)', ClosestMatchPosition2Hist(2,:)', 40, 'y');
end

figure(51); clf; title('Closest point Match Point'); hold on;
RawReScaledPosi2 = All_Ray(4:6,:).*repmat( ImgInfo(2).Model.Depth.RawDepth(:)'*DepthScale2, 3,1);
RawReScaledPosi1 = R*(All_Ray(1:3,:).*repmat( ImgInfo(1).Model.Depth.RawDepth(:)'*DepthScale1, 3,1)) + repmat(T, 1, length(All_Ray));
RawReScaledPosi2(:,IND2) = [];
RawReScaledPosi1(:,IND1) = [];
scatter3(RawReScaledPosi2(1,:)', RawReScaledPosi2(3,:)', RawReScaledPosi2(2,:)', 1*ones(1,size( RawReScaledPosi2,2)));
scatter3(RawReScaledPosi1(1,:)', RawReScaledPosi1(3,:)', RawReScaledPosi1(2,:)', 0.5*ones(1,size( RawReScaledPosi1,2)));
scatter3(ClosestMatchPosition2(1,:)', ClosestMatchPosition2(3,:)', ClosestMatchPosition2(2,:)', 40, 'g');
scatter3(ClosestMatchPosition2(1,:)', ClosestMatchPosition2(3,:)', ClosestMatchPosition2(2,:)', 40, 'b');
line( [ ClosestMatchPosition2(1,:); ClosestMatchPosition1(1,:)], ...
      [ ClosestMatchPosition2(3,:); ClosestMatchPosition1(3,:)], ...
      [ ClosestMatchPosition2(2,:); ClosestMatchPosition1(2,:)]);
line( [ MonoStichPosition2(1,:); MonoStichPosition1(1,:)], ...
      [ MonoStichPosition2(3,:); MonoStichPosition1(3,:)], ...
      [ MonoStichPosition2(2,:); MonoStichPosition1(2,:)]);
end  
return;

