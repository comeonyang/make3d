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
function []=SaveModelSatus( defaultPara, Name, ImgName, GroundLevel, R, T, Scale, FlagFirstPair)

% Function storage the GroundLevel, R, T, Scale of the ImgName in Name Model (Multi_Name.wrl)
[status dump] = system(['ls ' defaultPara.Fdir '/data/ModelStatus/' Name '.mat' ]);
if ~status && (~FlagFirstPair)
	load([ defaultPara.Fdir '/data/ModelStatus/' Name '.mat']);
else
    ModelStatus(1).ImgName = ImgName{1};
    ModelStatus(1).GroundLevel = GroundLevel;
    ModelStatus(1).R = defaultPara.R;
    ModelStatus(1).T = defaultPara.T;
    ModelStatus(1).Scale = defaultPara.Scale;
end

Ptr = length(ModelStatus);
ModelStatus(Ptr+1).ImgName = ImgName{2};
ModelStatus(Ptr+1).GroundLevel = GroundLevel;
ModelStatus(Ptr+1).R = R;
ModelStatus(Ptr+1).T = T;
ModelStatus(Ptr+1).Scale = Scale;

save([defaultPara.Fdir '/data/ModelStatus/' Name '.mat'],'ModelStatus');
return;
