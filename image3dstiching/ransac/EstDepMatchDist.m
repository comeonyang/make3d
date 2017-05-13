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
function [dist, inlierThreDist] = EstReProjError(X, R, T, D1, D2, lamda1, lamda2, disp)

% This function calculate the Estimated Reprojction Error for Match pairs
% According to heuristic
% Bigger Estimated Reprojction Error is lower the match is correct
% Input 
% 	X - calibrated corrdinate (normalize by depth)
% 	R - rotation matirx
% 	T - translation matrix (unit length)
%   D1/2 - depth information
%   lamda1/2 - triangulated depths
% Return
% 	dist - distribution define according to the heuristic
%	inlierThreDist - when EstDepMatchDist <= threDist it is a inlier

NumMatches = size(X,2);
Thre = Inf;
%  X1 = X(1:3,:);
%  X2 = X(4:6,:);
%  X2_2 = X2.*repmat(D2, 3, 1); 
%  X1_2 = R*X1.*repmat(D1, 3, 1); 

 ops = sdpsettings('solver','sedumi','verbose',1);
 a = sdpvar(1,1);
 b = sdpvar(1,1);
 F = set(a>=0)+set(b>=0);
%  sol = solvesdp(F,norm(a*X1_2(:) + repmat(T, NumMatches, 1)-
%  a*X2_2(:),2),ops);
 sol = solvesdp(F,norm( lamda1*a - D1, 1)+norm( lamda2*b - D2, 1),ops);
 a = double(a);
 b = double(b);
 EstReProjError = abs(D1./lamda2/a-lamda1./lamda2)+abs(D2./lamda1/b-lamda2./lamda1);
 inlierThreDist = find(EstReProjError <= Thre);
 if disp
    figure(6); 
    hist(EstReProjError(inlierThreDist),1000);
 end
 
 % fit exp distibution
 parmhat = expfit(EstReProjError(inlierThreDist));
 dist = exppdf(EstReProjError(inlierThreDist),parmhat);
 return;
