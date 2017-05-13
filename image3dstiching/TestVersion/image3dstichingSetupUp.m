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
% This function setup all the nessesary initialization for all the process in 
% image3dstiching folder or sub-folder

% Set Path
addpath(genpath('./'));

% Old function (should not be used anymore)
%addpath(genpath('../multipleImages')); % Old function

% general function
addpath(genpath('../third_party/lightspeed')); % lightspeed
addpath(genpath('../third_party/torr/')); % torr/
addpath(genpath('../third_party/vlutil')); % vlutil
addpath(genpath('../third_party/zisserman')); % zisserman

% Projective Reconstruction function
addpath(genpath('../third_party/missing-dat')); %missing-data

% feature matching function
addpath(genpath('../third_party/SURF-V1.0.8')); % SURF
addpath(genpath('../third_party/sift')); % SIFT
addpath(genpath('../third_party/kovesi')); % Ransac

% sfm function
addpath(genpath('../third_party/sba-1.3')); % Non-linear solver (sparse Levenberg-Marquardt algorithm)

% opt function
addpath(genpath('../third_party/opt')); % OPT solver

% LearningCode functions
addpath(genpath('../LearningCode')); % learned depth and plane parameter functions
% other need for LearningCode
addpath(genpath('../third_party/EdgeLinkLineSegFit/')); % EdgeLinkLineSegFit/

% initialize default Parameters
% foun from http://en.wikipedia.org/wiki/Earth_radius
defaultPara.ellip_equatorial_radius = 6378.135e3; % in meters
defaultPara.ellip_polar_radius = 6356.750e3;% in meters
defaultPara.VertVar = 0.25;% in percentage 

% Default camera intrinsic parameters (image pixel coordinate) from kyles' data
% assuming camera is horizontal pose (1704 by 2272)
defaultPara.fy = 2400.2091651084;
defaultPara.fx = 2407.3312729885838;
defaultPara.Ox = 1110.7122391785729;%1704/2;
defaultPara.Oy = 833.72104535435108;%2272/2;
defaultPara.InrinsicK1 = [ [defaultPara.fx 0 defaultPara.Ox];...
                   [0 -defaultPara.fy defaultPara.Oy];...
                   [0 0 1]];
defaultPara.InrinsicK2 = defaultPara.InrinsicK1;

