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
function [Model] = CleanUpMonoInfo2MultiInfo(defaultPara, ImgName)

% This function clean the loaded data and assign to specific structure
% Data loaded:
%	1) FitDepth 2) MedSup 3) MultiScaleSupTable 4) Ray 5) RayOri 6) Sup 
% 	7) SupNeighborTable 8) SupOri 9) depthMap 10) maskG 11) maskSky

% Return:
%	1) model.Depth - FitDepth
%		       - RawDepth
%	2) model.Ray
%	3) model.Sup
%	4) model.SupNeighborTable
%	5) model.MultiScaleSupTable

%add	Constrain.RayMatche Constrain.Depth_modified Constrain.SupMatched	
load( [defaultPara.ScratchFolder ImgName '/' ImgName '__AInfnew.mat']);

% assign the model value
Model.Depth.FitDepth = FitDepth;
Model.Depth.RawDepth = depthMap;
Model.Ray = Ray;
Model.Sup = Sup;
Model.maskG = maskG;
Model.maskSky = maskSky;
Model.SupNeighborTable = SupNeighborTable;
Model.MultiScaleSupTable = MultiScaleSupTable;

% add empty space holder
Model.Constrain.RayMatche = [];
Model.Constrain.Depth_modified = [];
Model.Constrain.SupMatched = [];

% clear thee rest
clear MedSup RayOri SupOri maskG maskSky;
return;

