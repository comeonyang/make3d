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
function [model]=AddMatch2Model(defaultPara, Wrlname, lamda, Xim, ImgInfo, ImgScale, ImgIndex1, ImgIndex2, PostFixStr, Error)

% This function load the previous mat file to get the model variables
% Then add new match constrain in the model.ConstrainOccluMatch

%load([defaultPar.Fdir '/data/' ]) % no need to load again ImgInfo is the latest version

% ============================================================================
% need Sup Ray Depth
% 1) Sup
NumMatches = size(Xim,2);
ScaleDepth = size(ImgInfo.Model.Depth.FitDepth);
[IND1] = ProjPosi2Mask( ImgScale, ScaleDepth, Xim);
SupMatched = ImgInfo.Model.Sup(IND1)'; % Sup
mask = SupMatched == 0; % erase sky points
SupMatched(mask)=[];

% 2) Ray
x_calib = inv(defaultPara.InrinsicK1)*[ Xim; ones(1, NumMatches)];
RayNormfactor = sqrt( sum(x_calib.^2,1));
RayMatched = (x_calib./(repmat( RayNormfactor,3,1)))'; % Ray
RayMatched(mask,:) = []; % erase sky points

% depth
% Notice: lamda is in the scale of local model already
Depth_modified = lamda.*RayNormfactor;
Depth_modified(:,mask) = []; % erase sky points
% =========================================================================

% Compensate the number of the OccluMatches to the OriMatches =============
% count No of matches point to determine the weights
NumMatches = size(ImgInfo.Model.Constrain.RayMatche,1);
NumOccluMatches = size(RayMatched,1);
% Ratio = ceil(NumMatches/NumOccluMatches);
Gap = NumMatches - NumOccluMatches;

if Gap > 0 
    if NumOccluMatches ~=1
        distMap = ( sum((repmat(Xim, [1 1 NumOccluMatches]) - repmat( permute(Xim,[ 1 3 2]), [1 NumOccluMatches 1])).^2,1));% do not do the sqrt to put in some nonliearity
        distMap(distMap ==0) = Inf;
        dist = min( distMap,[],3);
        distAffectByError = dist./Error;
        TotalInvDistAffectByError = sum(1./distAffectByError);
        NumDuplicate = ceil( Gap*(1./distAffectByError./TotalInvDistAffectByError));
%         NumDuplicate(1:NumOccluMatches)  = 1;
        for i = 1:NumOccluMatches
            RayMatched = [RayMatched; repmat(RayMatched(i,:), NumDuplicate(i),1)];
            Depth_modified = [ Depth_modified repmat(Depth_modified(:,i), 1, NumDuplicate(i))];
            SupMatched = [SupMatched; repmat(SupMatched(i,:), NumDuplicate(i),1)];        
        end   
    else
        NumDuplicate = Gap;
        RayMatched = [RayMatched; repmat(RayMatched, NumDuplicate,1)];
        Depth_modified = [ Depth_modified repmat(Depth_modified, 1, NumDuplicate)];
        SupMatched = [SupMatched; repmat(SupMatched, NumDuplicate,1)];
    end    
end    
% =======================================================================

% add new match constrain
ImgInfo.Model.ConstrainOccluMatch.RayMatche{ImgIndex2} = RayMatched;
ImgInfo.Model.ConstrainOccluMatch.Depth_modified{ImgIndex2} = Depth_modified;
ImgInfo.Model.ConstrainOccluMatch.SupMatched{ImgIndex2} = SupMatched;

% new model variable
model = ImgInfo.Model;
if defaultPara.Flag.FlagRefinementDisp
	disp('Storaging of Occlusion Match Correcting Single Modle in model.ConstrainOccluMatch');
end
save([defaultPara.Fdir '/data/' ImgInfo.ExifInfo.IDName '/' Wrlname '_' ImgInfo.ExifInfo.IDName '_' PostFixStr '.mat'],'model');

return;
