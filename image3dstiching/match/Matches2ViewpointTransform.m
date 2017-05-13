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
function []=Matches2ViewpointTransform(defaultPara, x1, x2, theta_hat, psi_hat, FlagDisp)

% This function generate the possible viewpoints of each matches
% after collecting all the possible viewpoints for all the matches
% we can count the most frequent viewpoints as the estimated viewpoints

% hope this approach will be robust to outlier and 
% the initial search Measured viewpoint theta_hat, psi_hat
% will make this approach better than (ransac then EstPose)

% restriction:
% we assume 1) only one axis of rotation
% 	    2) there is no height translation

% Input:
%	1) measured translation angel [cos(theta_hat) 0 sin(theta_hat)] = [x y z] ; y =0;
%	2) measured rotation angle [ [ cos(psi_hat)  0  sin(psi_hat) ];...
%				     [            0  1            0  ];...
%				     [ -sin(psi_hat) 0  cos(psi_hat) ]];

% Output:
%	ViewPoint :
%		   [ psi .......;...
%		     theta .....;...]
%                  the third axis is number of matches

% default paramter
psi_range = 30/180*pi;
psi_step = 0.1/180*pi;
NumMatches = size(x1,2);

% initialize variable
psi = (psi_hat -psi_range):psi_step:(psi_hat + psi_range);
NumSample = size(psi,2);

% calculate corresponding theta for each psi
theta = GenThetaFromPsi(x1, x2, psi);

%ViewPoint = cat(3, psi, theta);

% plot the 2D sample intensitiy figure
if FlagDisp
	figure;
	hold on;
	scatter(reshape(repmat(psi,NumMatches,1), [], 1), theta(:), 3);
end
