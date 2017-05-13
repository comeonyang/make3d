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
function [matches1 matches2] = GenMatches(defaultPara, ImgInfo, FlagDisp)

% This function generate Matches uing IMU and GPS info and Ransac and BA

% 1. Mono calulation or load the pre-calculated data ------------------------
ImgInfo(1).appendOpt = 0;
ImgInfo(2).appendOpt = 0;
[ ImgInfo] = SingleModelInfo(defaultPara, ImgInfo);

% initialize variables
Img1 = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
Img2 = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
[f1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % original features 
[f2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % original features 
[D1 IND] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1, ImgInfo(1).Model.Depth.FitDepth);
[D2 IND] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2, ImgInfo(1).Model.Depth.FitDepth);

% 1. extract Measuesd Position and orientation from GPS or IMU info
[MeasR MeasT] = InitPoseMeas(defaultPara, ImgInfo(1), ImgInfo(2));
 
% 2. Using Measures R and T and Mono-Depth to define mach search space constrain
% read in all surf features
[ Rc1, Rc2, ConS1, ConS2, ConSRough1, ConSRough2] = CalMatchSearchRegin(defaultPara, MeasR, MeasT, I1, I2, f1, f2, D1, D2, 1, FlagDisp);
Vector2Ipoint([Rc1; ConS1],[defaultPara.Fdir '/surf/'],['RConS_' Img1]);
Vector2Ipoint([Rc2; ConS2],[defaultPara.Fdir '/surf/'],['RConS_' Img2]);
Vector2Ipoint([ConSRough1],[defaultPara.Fdir '/surf/'],['RConSRough_' Img1]);
Vector2Ipoint([ConSRough2],[defaultPara.Fdir '/surf/'],['RConSRough_' Img2]);

% 3. Do match search with all combinations satisfying Constrain from 2) using ralative threshould
tic;
cd match
pwd
% system(['./surfMatchRConS.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' _ 0.3 0.7']);  
system(['./surfMatchRConS.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' Dense ' '0.3 0.7']);    % Parameter still need to be changed//Min
cd ..
toc

% 4. Ransac
[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, [ defaultPara.Type 'Dense'], 1, 1);
if isempty(matches)
   disp('Zeros matches');
   matches1 = matches(1,:);
   matches2 = matches(2,:);   
   return;
end    
[D1 IND1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1(:,matches(1,:)), ImgInfo(1).Model.Depth.FitDepth);
[D2 IND2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2(:,matches(2,:)), ImgInfo(1).Model.Depth.FitDepth);
%figure(11);  plotmatches(I1,I2,f1, f2,matches, 'Stacking','v','Interactive', FlagDisp); title('SurfMatch')
%saveas(11,[defaultPara.ScratchFolder Img1 '_' Img2 'SimpleSurfMatch'],'jpg');
[F, inliers, NewDist, fail]=GeneralRansac(defaultPara, f1, f2, matches, D1, D2);
figure(12);  plotmatches(I1,I2,f1, f2,matches(:,inliers), 'Stacking', 'v', 'Interactive', FlagDisp);
saveas(12,[defaultPara.ScratchFolder Img1 '_' Img2 'AfterRansac'],'jpg');
close 12;

% *** Stop maunally to pick out the bad matches*** -----------------
matches = matches(:,inliers);
if isempty(matches)
   disp('Zeros matches');
   matches1 = matches(1,:);
   matches2 = matches(2,:);
   return;
end  

% x_calib = [ inv(defaultPara.InrinsicK1)*[ f1(:,matches(1,:)); ones(1,length(matches))];...
% 	      inv(defaultPara.InrinsicK2)*[ f2(:,matches(2,:)); ones(1,length(matches))]]; 
% [ lamda1 lamda2] = triangulation( defaultPara, MeasR(1:3,:), MeasT(1:3), x_calib);
% %    end
% X_obj_1 = x_calib(1:3,:).*repmat(lamda1, 3, 1);
% X_obj_2 = MeasR(4:6,:)*(x_calib(4:6,:).*repmat(lamda2, 3, 1)) + repmat(MeasT(4:6), 1, length(matches));
% X_obj = (X_obj_1+X_obj_2)/2;
% %end

% 5. Bundle Adjustment
% [R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, MeasR, MeasT, [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj, ImgInfo, 1);
% outlier_thre1 = prctile(dist1_BA,90);
% outlier_thre2 = prctile(dist2_BA,90);
% Outlier = dist1_BA > outlier_thre1 | dist2_BA > outlier_thre2;
% lamda1(Outlier) = [];
% lamda2(Outlier) = [];
% X_obj_BA(:,Outlier) = [];
% x_calib(:,Outlier) = [];
% matches(:, Outlier) = [];
% % [R T X_obj_BA X_im_BA dist1_BA dist2_BA]=SparseBAWraper(defaultPara, R, T, [f1(:,matches(1,:)); f2(:,matches(2,:))], X_obj_BA, ImgInfo, 1);
% figure(13); plotmatches(I1,I2,f1, f2,matches, 'Stacking', 'v', 'Interactive', FlagDisp);title('after BA clean once');
% saveas(13,[defaultPara.ScratchFolder Img1 '_' Img2 'AfterBA'],'jpg');

matches1 = matches(1,:);
matches2 = matches(2,:);
return;
