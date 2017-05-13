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
function []=Indep_Mono_Model_CheckRT(defaultPara, ImgInfo, R, T,Scale)
% []=Indep_Mono_Model_CheckRT(defaultPara, ImgInfo, R, T)

Default.fy = 2400.2091651084;
Default.fx = 2407.3312729885838;
Default.Ox = 1110.7122391785729;%2272/2; %
Default.Oy = 833.72104535435108;%1704/2; %
Default.a_default = 2272/Default.fx;
Default.b_default = 1704/Default.fy; %0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
Default.Ox_default = 1-Default.Ox/2272;%0.489272914; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;

% This function render the 2 independent mono-model jointly using RT info
Img1 = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
Img2 = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift

[WrlPosition1 PositionTex1] = PreWrlData(defaultPara.InrinsicK1, ImgInfo(1).Model.Depth.FitDepth*Scale(1), I1, R.a, T.a);
WrlFacestHroiReduce( WrlPosition1, double(PositionTex1), ImgInfo(1).Model.Sup, Img1, ['_MonoJointlyBA'], defaultPara.OutPutFolder, 0, 0);
                
system(['cp ' defaultPara.Fdir '/jpg/' Img1 '.jpg '  defaultPara.OutPutFolder Img1 '.jpg']);
[WrlPosition2 PositionTex2] = PreWrlData(defaultPara.InrinsicK2, ImgInfo(2).Model.Depth.FitDepth*Scale(2), I2, R.b, T.b);
WrlFacestHroiReduce(WrlPosition2, PositionTex2, ImgInfo(2).Model.Sup, Img2, [ '_MonoJointlyBA'], defaultPara.OutPutFolder, 0, 0);
system(['cp ' defaultPara.Fdir '/jpg/' Img2 '.jpg ' defaultPara.OutPutFolder Img2 '.jpg']);

return;
