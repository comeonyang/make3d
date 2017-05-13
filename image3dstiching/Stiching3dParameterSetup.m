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
% This script initialize all the parameter
% category include:
% 	1) GPS parameter
%	2) camera parameter
%	3) mono-info parameter
%	4) Specified folder
% 	5) Features and Matching Parameter
% 	6) Optimization parameter

% 1. GPS parameter
% foun from http://en.wikipedia.org/wiki/Earth_radius
defaultPara.ellip_equatorial_radius = 6378.137e3;%6378.135e3; % in meters
defaultPara.ellip_polar_radius = 6356.75231e3;%6356.750e3;% in meters

% 2. camera parameter
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
defaultPara.R = eye(3);
defaultPara.T = zeros(3,1);
defaultPara.Scale = 1;

% 3. mono-info parameter
defaultPara.ParaFolder = '/afs/cs/group/reconstruction3d/scratch/Para/';
defaultPara.taskName = '';
defaultPara.Flag = Flag;
defaultPara.Flag.DisplayFlag = 0;
defaultPara.Flag.IntermediateStorage = 0;
defaultPara.Flag.FeaturesOnly = 0;
defaultPara.Flag.NormalizeFlag = 1;
defaultPara.Flag.BeforeInferenceStorage = 1;
defaultPara.Flag.NonInference = 0;
defaultPara.Flag.AfterInferenceStorage = 1;
defaultPara.Flag.FeatureStorage = 0;



% 4. Specified folder
defaultPara.Fdir = Fdir;
defaultPara.ScratchFolder = [Fdir '/data/'];
%defaultPara.ScratchFolder = '/afs/cs/group/reconstruction3d/scratch/temp';
%defaultPara.OutPutFolder = '/afs/cs/group/reconstruction3d/scratch/3DmodelMultipleImage/';
defaultPara.OutPutFolder = [Fdir '/Wrl/'];%defaultPara.ScratchFolder;

% 5. Features and Matching Parameter
defaultPara.Type = '_RConS';

% 6. Optimization parameter
defaultPara.TriLeastSquare = 1;
defaultPara.opt = sdpsettings('solver','sedumi','cachesolvers',1,'removeequalities',2,'verbose',0);

% 7. Rendering Parameter
defaultPara.Wrlname = 'Multi';
defaultPara.RenderFlag = true;

% surf and matching parameter
defaultPara.AbsThre = 0.2;
defaultPara.RatioThre = 0.6;
defaultPara.VertVar = 0.25;% in percentage 
defaultPara.MaxRatio = 300; %Min changed
defaultPara.MinimumNumMatches = 50;
defaultPara.ReProjErrorThre = 1;
defaultPara.NegativeDepthTolerence = 0.5;
defaultPara.surf.SparseThre = 30000;
defaultPara.surf.DenseThre = 15000;
defaultPara.OccluDistThre = 10; % occlusion bigger than defaultPara.OccluDistThre then triangulate
defaultPara.radius2imageSizeRatio = 100; % used by CalMatchDensityWeights.m

% Ransac
defaultPara.MAXEnsembleSamples = 10000;

% Model parameter
defaultPara.InitialDepth = 'RawDepth';
defaultPara.FarestDist = 1e3; % Do not model structure lkm away from the camera
defaultPara.Closestdist = 1; % Do not model structure 1m away from the camera 

% other
%defaultPara.FlagFirstPair = 1;
defaultPara.LastImgFlag = 0;

% Refinement Thre
defaultPara.CoeffMThre = 0.8;
defaultPara.ResidualThre = 5;
defaultPara.coeffratioThre = 0.8;
defaultPara.PostFixStrAfter = 'PairNew';
defaultPara.MaxUniqueRatio = 50;
defaultPara.OcclusionRatioThre = 1.1;
