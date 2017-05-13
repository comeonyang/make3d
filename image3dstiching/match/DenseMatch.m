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
function [Pair]=DenseMatch(defaultPara, R, T, ImgInfo)

% This function search for denser mach given pretty accurate R T

Img1 = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
Img2 = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
[f1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % original features 
[f2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % original features 
[D1] = PorjPosi2Depth(size(I1), size(ImgInfo(1).Model.Depth.FitDepth), f1, ImgInfo(1).Model.Depth.FitDepth);
[D2] = PorjPosi2Depth(size(I2), size(ImgInfo(2).Model.Depth.FitDepth), f2, ImgInfo(1).Model.Depth.FitDepth);

% 1. Using BA's R and T and Scaled Mono-Depth to define match search space constrain
% read in all surf features
defaultPara.VertVar = 0.02;
defaultPara.MaxRatio = 5;
[ Rc1, Rc2, ConS1, ConS2, ConSRough1, ConSRough2] = CalMatchSearchRegin(defaultPara, [R; R'], [T; -R'*T], I1, I2, f1, f2, D1, D2, 1, 0);
Vector2Ipoint([Rc1; ConS1],[defaultPara.Fdir '/surf/'],['RConS_' Img1]);
Vector2Ipoint([Rc2; ConS2],[defaultPara.Fdir '/surf/'],['RConS_' Img2]);
Vector2Ipoint([ConSRough1],[defaultPara.Fdir '/surf/'],['RConSRough_' Img1]);
Vector2Ipoint([ConSRough2],[defaultPara.Fdir '/surf/'],['RConSRough_' Img2]);

% 2. Do match search with all combinations satisfying Constrain from 2) using ralative threshould
tic
cd match
system(['./surfMatchRConS.sh ' defaultPara.Fdir ' ' Img1 ' ' Img2 ' Dense ' '0.1 0.3']);    
cd ..
toc

[f1, f2, matches] = readSurfMatches(Img1, Img2, defaultPara.Fdir, [ defaultPara.Type 'Dense'], 1, 1);
figure; plotmatches(I1,I2,f1, f2,matches, 'Stacking', 'v', 'Interactive', 2);
f1 = f1(:,matches(1,:));
f2 = f2(:,matches(2,:));

% % 3. triangulation
% x_calib = [ inv(defaultPara.InrinsicK1)*[ f1; ones(1,length(f1))];...
% 	      inv(defaultPara.InrinsicK2)*[ f2; ones(1,length(f2))]];
% [ Pair.depth1 Pair.depth2] = triangulation( defaultPara, R, T, x_calib);
% X_obj_1 = x_calib(1:3,:).*repmat(Pair.depth1, 3, 1);
% X_obj_2 = R'*(x_calib(4:6,:).*repmat(Pair.depth2, 3, 1)) + repmat(-R'*T, 1, length(f1));
% Structure.X_obj = (X_obj_1+X_obj_2)/2;
Pair.Xim = [f1; f2];
Pair.R = R;
Pair.T = T;
return;
