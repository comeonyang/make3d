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
function [inlier, Residual] = EpipoPrune(defaultPara, Pair, Xim, ScaleImg)

% this function Prune out the matches that are far from the epipolae line
% initalize parameter
EpipolarThre = 0.05; % test 0.005 in Quid case still some error

% Building F
T1_hat = [[0 -Pair.T(3) Pair.T(2)];...
            [Pair.T(3) 0 -Pair.T(1)];...
            [-Pair.T(2) Pair.T(1) 0]];
F1 = inv(defaultPara.InrinsicK2)'*T1_hat*Pair.R*inv(defaultPara.InrinsicK1);

% pick Xim(1:2,:) check Xim(3:4,:)
EpipolarLinePara = F1*[ Xim(1:2,:); ones(1, size(Xim(1:2,:),2))];
EpipolarLinePara = EpipolarLinePara ./ repmat( sqrt( sum(EpipolarLinePara(1:2,:).^2,1)),3,1);
Mask1 = abs( sum( EpipolarLinePara .*  [ Xim(3:4,:); ones(1, size(Xim(3:4,:),2))])) <= EpipolarThre*max(ScaleImg);

% pick Xim(3:4,:) check Xim(1:2,:)
EpipolarLinePara = F1'*[ Xim(3:4,:); ones(1, size(Xim(1:2,:),2))];
EpipolarLinePara = EpipolarLinePara ./ repmat( sqrt( sum(EpipolarLinePara(1:2,:).^2,1)),3,1);
Residual = abs( sum( EpipolarLinePara .*  [ Xim(1:2,:); ones(1, size(Xim(1:2,:),2))]));
Mask2 = Residual <= EpipolarThre*max(ScaleImg);

inlier = Mask1 & Mask2;

return;
