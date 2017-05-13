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
            [Matches1 CoeffM1 Inliers1]=CorrolationMatch( defaultPara, Pair, I1, I2, PointPix1(:,Mask1), POriReprojM1(:,Mask1), FieldOccluPix1(:,Mask1), [minRatio maxRatio]);
            save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter 'FirstCorrMatches.mat'],'Matches1','CoeffM1','Inliers1');
		    Pair2_1.R = Pair.R';
	   	    Pair2_1.T = -Pair.R'*Pair.T;
                    [Matches2 CoeffM2 Inliers2]=CorrolationMatch( defaultPara, Pair2_1, I2, I1, PointPix2(:,Mask2), POriReprojM2(:,Mask2), FieldOccluPix2(:,Mask2), [1/maxRatio 1/minRatio]);
            save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter 'SecondCorrMatches.mat'],'Matches2','CoeffM2','Inliers2');

if true
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

  		    % use Coeff as threshould to filter out error matches
		    CoeffMask = CoeffM(1,:) > defaultPara.CoeffMThre;
		    [inlier, Residual] = EpipoPrune(defaultPara, Pair, Matches, ImgScale1);
		    EpipolarResidualMask = Residual < defaultPara.ResidualThre;
		    CoeffRationMask = CoeffM(2,:)./CoeffM(1,:) < defaultPara.coeffratioThre;
		    Mark = CoeffMask & CoeffRationMask & EpipolarResidualMask;
		    % =======================================================================
                    if ~isempty(Matches)
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
                    save([defaultPara.Fdir '/data/' Img1 '_' Img2 '_' PostFixStrAfter 'TestScript.mat'],'Matches','CoeffM','Error','Mark');
                    end

                    % Storage the New Matches                   
		    if defaultPara.Flag.FlagRefinementDisp
			disp('Storaging Occlusion Surf Features Matches');
		    end
        
end            
