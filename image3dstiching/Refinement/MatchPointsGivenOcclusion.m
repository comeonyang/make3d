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
function [matches fail] = MatchPointsGivenOcclusion( defaultPara, ImgInfo1, ImgInfo2, ImgScale1, ImgScale2, Img1, Img2, Img1Index, Img2Index, ...
                    Pair, GlobalScale, Wrlname, PostFixStrAfter,...
					POriReprojM1, FieldOccluPix1, FaceSetPickedIND1, ...
					POriReprojM2, FieldOccluPix2, FaceSetPickedIND2, FlagEarlyStopMatchPointsGivenOcclusion)
% initialize parameters
displayFlag = 0;
fail = 0;
RefineCorrSpace = 10;
depthratioMin = 0.01;
depthratioMax = 100;

if nargin < 20
	FlagEarlyStopMatchPointsGivenOcclusion = 0;
end
% 1)Prepare the Constrain to run the matching 

% load surf Features
[f1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % original features
[f2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % original features

% initialize the Rc ConS ConSRough
NumSurF1 = length(f1);
NumSurF2 = length(f2);
Rc1 = [ ones(1,NumSurF1); zeros(1,NumSurF1); ones(1,NumSurF1); zeros(1,NumSurF1)];
Rc2 = [ ones(1,NumSurF2); zeros(1,NumSurF2); ones(1,NumSurF2); zeros(1,NumSurF2)];
ConS1 = zeros(4,NumSurF1);
ConS2 = zeros(4,NumSurF2);
ConSRough1 = ConS1;
ConSRough2 = ConS2;
ConS1_4points = zeros(8,NumSurF1);
ConS2_4points = zeros(8,NumSurF2);
AllPOriReprojM1 = zeros(2,NumSurF1);
AllPOriReprojM2 = zeros(2,NumSurF2);
AllFieldOccluPix1 = zeros(2,NumSurF1);
AllFieldOccluPix2 = zeros(2,NumSurF2);

% calculate constrain for effective points
if ~isempty(FaceSetPickedIND1)
    [ Rc1(:,FaceSetPickedIND1) ConS1(:,FaceSetPickedIND1) ConSRough1(:,FaceSetPickedIND1) ConS1_4points(:,FaceSetPickedIND1)] = ...
        EndPoint2BoxConS(defaultPara, ImgScale1(1), ImgScale1(2), POriReprojM1, FieldOccluPix1, 1);
	AllPOriReprojM1(:,FaceSetPickedIND1) = POriReprojM1;
	AllFieldOccluPix1(:,FaceSetPickedIND1) = FieldOccluPix1;
end
if ~isempty(FaceSetPickedIND2)
    [ Rc2(:,FaceSetPickedIND2) ConS2(:,FaceSetPickedIND2) ConSRough2(:,FaceSetPickedIND2) ConS2_4points(:,FaceSetPickedIND2)] = ...
        EndPoint2BoxConS(defaultPara, ImgScale2(1), ImgScale2(2), POriReprojM2, FieldOccluPix2, 1);
	AllPOriReprojM2(:,FaceSetPickedIND2) = POriReprojM2;
	AllFieldOccluPix2(:,FaceSetPickedIND2) = FieldOccluPix2;
end

T1_hat = [[0 -Pair.T(3) Pair.T(2)];...
    [Pair.T(3) 0 -Pair.T(1)];...
    [-Pair.T(2) Pair.T(1) 0]];
F = inv(defaultPara.InrinsicK2)'*T1_hat*Pair.R*inv(defaultPara.InrinsicK1);
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
if displayFlag

	figure;
	dispMatchSearchRegin(I1, I2, [f1; ones(1,NumSurF1)], [f2; ones(1,NumSurF2)], ConS1_4points, ConS2_4points, F, ...
		AllPOriReprojM1, ones(1,NumSurF1), AllFieldOccluPix1, ones(1,NumSurF1), ...
		AllPOriReprojM2, ones(1,NumSurF2), AllFieldOccluPix2, ones(1,NumSurF2), ...
		1, 'Stacking', 'h', 'Interactive', 0);
end

if ~( isempty(FaceSetPickedIND1)&&isempty(FaceSetPickedIND2)) % only if not both FaceSetPickedIND is empty then find the matches
	% write the constrain into data
	Vector2Ipoint([Rc1; ConS1],[defaultPara.Fdir '/surf/'],['RConS_' Img1]);
	Vector2Ipoint([Rc2; ConS2],[defaultPara.Fdir '/surf/'],['RConS_' Img2]);
	Vector2Ipoint([ConSRough1],[defaultPara.Fdir '/surf/'],['RConSRough_' Img1]);
	Vector2Ipoint([ConSRough2],[defaultPara.Fdir '/surf/'],['RConSRough_' Img2]);

        %======================= debug only
%        save([defaultPara.Fdir '/data/PreOcclusionDetect.mat'],'Rc1','ConS1','Rc2','ConS2','ConSRough1','ConSRough2');
%        return;
	%==================================

	% run time consuming matching code
	tic;
	cd match
		system(['./surfOccluMatch.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' OccluDense ' '0.2 0.6']);    % Parameter still need to be changed//Min
	cd ..
	toc

	% Readin matching result
	[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, [ defaultPara.Type 'OccluDense'], 1, 1, 3);
        % in 'OccluDense' cases matches is N by 3, the last column is Ratio
        if isempty( matches)
            disp( 'Zeros Surf Occlusion Matches');
            fail = 1;
            return; 
        end
        matches = matches(1:2,:);

else
	% no matches, means no need to storage new ConstrainOccluMatch
	matches = [];
end

if displayFlag
	figure(200); plotmatches(I1,I2,f1, f2,matches, 'Stacking','v','Interactive', 3);
	saveas(200,[ defaultPara.ScratchFolder Img1 '_' Img2 '_OccluMatches'],'jpg');
end

if FlagEarlyStopMatchPointsGivenOcclusion
	matches = [f1(:,matches(1,:)); f2(:,matches(2,:))];
    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'f1','f2','matches');
	return;
end

% 2) 	Process the matches

% Pruning by epipolarline
if ~isempty(matches)
	[inlier Residual] = EpipoPrune(defaultPara, Pair, [f1(:,matches(1,:)); f2(:,matches(2,:))], (ImgScale1+ImgScale2)/2);
	matches = matches(:,inlier);

	if defaultPara.Flag.FlagCorrRefinement
		% Fineer Search of the close by the SurfMatches Features by Corrolation Matches ========Min Added July 13th
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

		% construction epipolar line unit vector
		EpipolarUnitVector1 = F*[ f1(:,matches(1,:)); ones(1,size(matches,2))];	
		EpipolarUnitVector2 = F'*[ f2(:,matches(2,:)); ones(1,size(matches,2))];
		EpipolarUnitVector1 = EpipolarUnitVector1([2 1],:);	
		EpipolarUnitVector2 = EpipolarUnitVector2([2 1],:);	
		EpipolarUnitVector1 = EpipolarUnitVector1./repmat( sqrt( sum( EpipolarUnitVector1.^2, 1)), 2, 1);
		EpipolarUnitVector2 = EpipolarUnitVector2./repmat( sqrt( sum( EpipolarUnitVector2.^2, 1)), 2, 1);
	
		[Matches1 CoeffM1 Inliers1]=CorrolationMatch( defaultPara, Pair, I1, I2, f1(:,matches(1,:)), ...
			 f2(:,matches(2,:)) + EpipolarUnitVector1*RefineCorrSpace, ...
			 f2(:,matches(2,:)) - EpipolarUnitVector1*RefineCorrSpace, [minRatio maxRatio],[1 2]);
		Pair2_1.R = Pair.R';
		Pair2_1.T = -Pair.R*Pair.T;
		[Matches2 CoeffM2 Inliers2]=CorrolationMatch( defaultPara, Pair2_1, I2, I1, f2(:,matches(2,:)), ...
			 f1(:,matches(1,:)) + EpipolarUnitVector2*RefineCorrSpace, ...
			 f1(:,matches(1,:)) - EpipolarUnitVector2*RefineCorrSpace, [minRatio maxRatio],[1 2]);
		Matches1 = Matches1(:,Inliers1);
		Matches2 = Matches2(:,Inliers2);
		CoeffM1 = CoeffM1(Inliers1);
		CoeffM2 = CoeffM2(Inliers2);

		% Check if the Matches are not mutual discard the one with less Coeff(Cross-Corrolation value) ===============
		Matches = [ Matches1 [Matches2(3:4,:); Matches2(1:2,:)]];
		CoeffM = [ CoeffM1 CoeffM2];
		% Min used different algorithm than SurFeature Matches
		[Inliers] = CleanMatch(Matches, CoeffM); % choose the matches with higher Coeff is the matches is not mutual
		Matches = Matches(:,Inliers);
		CoeffM = CoeffM(:,Inliers);
		if defaultPara.Flag.FlagRefinementDisp
			figure; plotmatches(I1,I2,Matches(1:2,:), Matches(3:4,:),repmat(1:size(Matches,2),2,1), 'Stacking','v','Interactive', 3);
		end
		
		    % use Coeff as threshould to filter out error matches
		    Mask = CoeffM > defaultPara.CoeffMThre;
		    [inlier, Residual] = EpipoPrune(defaultPara, Pair, Matches, ImgScale1);
		    Mark = Mask & Residual < defaultPara.ResidualThre;
		    f1 = Matches(1:2,Mark);
		    f2 = Matches(3:4,Mark);
		    matches = repmat( 1:sum(Mark),2,1);
	end
	% ======================================================================================

	if displayFlag
		figure(201); plotmatches(I1,I2,f1, f2,matches, 'Stacking','v','Interactive', 3);
		saveas(201,[ defaultPara.ScratchFolder Img1 '_' Img2 '_OccluMatchesPrune'],'jpg');
	end
end

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
        AddMatch2Model(defaultPara, Wrlname, lamda1, f1(:,matches(1,:)), ImgInfo1, ImgScale1, Img1Index, Img2Index, PostFixStrAfter, Error);
    	AddMatch2Model(defaultPara, Wrlname, lamda2, f2(:,matches(2,:)), ImgInfo2, ImgScale2, Img2Index, Img1Index, PostFixStrAfter, Error);
end

% Storage the New Matches
if defaultPara.Flag.FlagRefinementDisp
	disp('Storaging Occlusion Surf Features Matches');	
end

save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter '.mat'],'f1','f2','matches', 'fail');

return;
