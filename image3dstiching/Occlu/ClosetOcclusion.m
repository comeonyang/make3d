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
function [Pocclu] = ClosetOcclusion(defaultPara, ScaleImg, Pair, ModelInfoTarget, TargetPointPix, INDTarget, ModelInfoField, FieldPointPix, INDField)
% function return the closest occlusion point on specific ray TargetPointPix

% ==========================================

% PlaneParaMeter Scaling and Transform to Target viewpoint
PlaneParaField = ModelInfoField.Model.PlaneParaInfo.PlanePara;
PlaneParaField = PlaneParaField( :,ModelInfoField.Model.PlaneParaInfo.Sup2Para(ModelInfoField.Model.PlaneParaInfo.SupEpand(INDField)));
PlaneParaField2Target = PlaneParaField'*Pair.R./repmat(1-PlaneParaField'*Pair.T,1,3); % (n by 3


% Find intersection with ray and the PlaneParaMeter
OccluDepthTarget = 1./(PlaneParaField2Target*ModelInfoTarget.Model.Ray(:,INDTarget))'; % (3 by n)
InterSectP = repmat(ModelInfoTarget.Model.Ray(:,INDTarget),1,length( OccluDepthTarget)).*repmat(OccluDepthTarget,3,1);

% Reproj the intersect point to Field viewpoint, check if in the superpixel that planes belong to
NumINDField = length( INDField);
ReProjField = defaultPara.InrinsicK1*(Pair.R*InterSectP + repmat( Pair.T, 1,NumINDField)); 
ReProjField = ReProjField(1:2,:)./repmat( ReProjField(3,:),2,1);
ScaleDepth = size(ModelInfoTarget.Model.Depth.FitDepth);
[ ReProjFieldIND] = ProjPosi2Mask( ScaleImg, ScaleDepth, ReProjField);
OccluMark = ModelInfoField.Model.Sup(ReProjFieldIND) == ModelInfoField.Model.Sup(INDField);


% find the closest occlued point
OccluDepthTarget = OccluDepthTarget(OccluMark);
[V I ]= min(OccluDepthTarget);
Pocclu=ReProjField(:,OccluMark);
Pocclu = Pocclu(:,I);

if V > ModelInfoTarget.Model.PlaneParaInfo.FitDepth(INDTarget)
    Pocclu = [];
end

return;


