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
function [ Matches]=CorrMatchPointsGivenOcclusion( ...
		defaultPara, Pair, GlobalScale,...
		POriReprojM1, FieldOccluPix1,PointPix1,...
		POriReprojM2, FieldOccluPix2,PointPix2,...
		ImgInfo1, ImgInfo2, ImgScale1, ImgScale2, ...
		FlagEarlyStopCorrMatchPointsGivenOcclusion)

% This function perform correlation matches given Occlusion info

if nargin < 12
	FlagEarlyStopCorrMatchPointsGivenOcclusion = 0;
end
PostFixStrAfter = 'CorrMatches';
Img1 = ImgInfo1.ExifInfo.IDName
Img2 = ImgInfo2.ExifInfo.IDName
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift		

depthratioMin = 0.01;
depthratioMax = 100;

    if ~isempty(Pair.lamda)
	depthratio = Pair.lamda(1,:)./Pair.lamda(2,:);
	% Min add to remove outliers (should be already removed when doing PoseEst.m)
	Inliers = depthratio > depthratioMin & depthratio < depthratioMax;
	depthratio = depthratio(Inliers);
	if ~isempty(depthratio)
		maxRatio = max(depthratio);
		minRatio = min(depthratio);
	else
		maxRatio = depthratioMax;
		minRatio = depthratioMin;
	end
    else
	maxRatio = max([GlobalScale(1)/GlobalScale(2) GlobalScale(2)/GlobalScale(1)]);
	minRatio = 1/maxRatio;
    end

    [Matches1 CoeffM1 Inliers1]=CorrolationMatch( defaultPara, Pair, I1, I2, PointPix1, POriReprojM1, FieldOccluPix1, [minRatio maxRatio]);
    if defaultPara.Flag.FlagRefinementDisp
	disp('Start Reverse CorrolationMatches Hold on');
% 	whos
    end

    Pair2_1.R = Pair.R';
    Pair2_1.T = -Pair.R'*Pair.T;
    [Matches2 CoeffM2 Inliers2]=CorrolationMatch( defaultPara, Pair2_1, I2, I1, PointPix2, POriReprojM2, FieldOccluPix2, [1/maxRatio 1/minRatio]);
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
	[InliersReverse] = CleanMatch(Matches([3 4 1 2],Inliers), CoeffM(1,Inliers)); % choose the matches with higher Coeff is the matches is not mutual
	Inliers = Inliers(InliersReverse);
	Matches = Matches(:,Inliers);
	CoeffM = CoeffM(:,Inliers);

    if defaultPara.Flag.FlagRefinementDisp
	figure; plotmatches(I1,I2,Matches(1:2,:), Matches(3:4,:),repmat(1:size(Matches,2),2,1), 'Stacking','v','Interactive', 3);
    end

%    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter 'BeforeFiltering.mat'],'Matches','CoeffM');

    % use Coeff as threshould to filter out error matches
    CoeffMask = CoeffM(1,:) > defaultPara.CoeffMThre;
    [inlier, Residual] = EpipoPrune(defaultPara, Pair, Matches, ImgScale1);
    EpipolarResidualMask = Residual < defaultPara.ResidualThre;
    CoeffRationMask = CoeffM(2,:)./CoeffM(1,:) < defaultPara.coeffratioThre;
    CoeffRationMask( CoeffM(1,:) == 0) = false;% CoeffM(1,:) == 0 means not satisfy the epipolar constrain
    Mark = CoeffMask & CoeffRationMask & EpipolarResidualMask;

%    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'Matches','CoeffM','Mark');

if FlagEarlyStopCorrMatchPointsGivenOcclusion
	Matches = Matches(:,Mark);
    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'Matches','CoeffM');
	return;
end
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
	AddMatch2Model(defaultPara, defaultPara.Wrlname, lamda1, Matches(1:2,Mark), ImgInfo1, ImgScale1, Img1Index, Img2Index, PostFixStrAfter, Error);
	AddMatch2Model(defaultPara, defaultPara.Wrlname, lamda2, Matches(3:4,Mark), ImgInfo2, ImgScale2, Img2Index, Img1Index, PostFixStrAfter, Error);
	% Storage the New Matches                   
%	if defaultPara.Flag.FlagRefinementDisp
		disp('Storaging Occlusion Surf Features Matches');
		save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'Matches','CoeffM','Error','Mark');
%	end
    end

return;
