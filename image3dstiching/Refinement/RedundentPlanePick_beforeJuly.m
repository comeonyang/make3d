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
function [model1 model2] = RedundentPlanePick(defaultPara, PointPix1, POriReprojM1, PointPix2, POriReprojM2, ImgInfo1, ImgInfo2, ScaleImg1, ScaleImg2)

IND1 =  find( sum(PointPix1,1) ~= 0);
IND2 =  find( sum(PointPix2,1) ~= 0);
ScaleDepth1 = size(ImgInfo1.Model.Depth.FitDepth);
ScaleDepth2 = size(ImgInfo2.Model.Depth.FitDepth);

% find the index of the sup given the point and occlusion position in image pixels
[Img1IndexTarget] = ProjPosi2Mask( ScaleImg1, ScaleDepth1, PointPix1(:,IND1));
[Img2IndexField] = ProjPosi2Mask( ScaleImg2, ScaleDepth2, POriReprojM1(:,IND1));
    % clean Target or Field map to sky (SupExpand() ==0)
    SkymaskTarget = ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexTarget) == 0;
    SkymaskField = ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexField) == 0;
    Img1IndexTarget(SkymaskTarget | SkymaskField) = [];
    Img2IndexField(SkymaskTarget | SkymaskField) = [];    
    
[Img2IndexTarget] = ProjPosi2Mask( ScaleImg2, ScaleDepth2, PointPix2(:,IND2));
[Img1IndexField] = ProjPosi2Mask( ScaleImg1, ScaleDepth1, POriReprojM2(:,IND2));
 % clean Target or Field map to sky (SupExpand() ==0)
    SkymaskTarget = ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexTarget) == 0;
    SkymaskField = ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexField) == 0;
    Img2IndexTarget(SkymaskTarget | SkymaskField) = [];
    Img1IndexField(SkymaskTarget | SkymaskField) = [];


% take out the Plane Parameter to check the orientation
PlanePara1Target = ImgInfo1.Model.PlaneParaInfo.PlanePara(:,ImgInfo1.Model.PlaneParaInfo.Sup2Para(...
                   ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexTarget)));
PlanePara1Target = PlanePara1Target./repmat(sqrt( sum( PlanePara1Target.^2,1)), 3, 1);
PlanePara1Field = ImgInfo2.Model.PlaneParaInfo.PlanePara(:,ImgInfo2.Model.PlaneParaInfo.Sup2Para(...
                   ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexField)));
PlanePara1Field = PlanePara1Field./repmat(sqrt( sum( PlanePara1Field.^2,1)), 3, 1);

PlanePara2Target = ImgInfo2.Model.PlaneParaInfo.PlanePara(:,ImgInfo2.Model.PlaneParaInfo.Sup2Para(...
                    ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexTarget)));
PlanePara2Target = PlanePara2Target./repmat(sqrt( sum( PlanePara2Target.^2,1)), 3, 1);
PlanePara2Field = ImgInfo1.Model.PlaneParaInfo.PlanePara(:,ImgInfo1.Model.PlaneParaInfo.Sup2Para(...
                    ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexField)));
PlanePara2Field = PlanePara2Field./repmat(sqrt( sum( PlanePara2Field.^2,1)), 3, 1);


% Pick out the Plane has a bigger normalized PlanPrameter of the Z component
Img1IndexTargetOutliers = abs(PlanePara1Target(3,:)) < abs(PlanePara1Field(3,:)); 
Img2IndexFieldOutliers = ~Img1IndexTargetOutliers;
Img2IndexTargetOutliers = abs(PlanePara2Target(3,:)) < abs(PlanePara2Field(3,:)); 
Img1IndexFieldOutliers = ~Img2IndexTargetOutliers;


% Modify the Sup to 0 so that it won't rendered
tempSup = ImgInfo1.Model.PlaneParaInfo.SupOri;
tempSup( Img1IndexTarget( Img1IndexTargetOutliers)) = 0;
tempSup(:,end) = ImgInfo1.Model.PlaneParaInfo.SupOri(:,end);
tempSup(end,:) = ImgInfo1.Model.PlaneParaInfo.SupOri(end,:);
ImgInfo1.Model.PlaneParaInfo.SupOri = tempSup;

tempSup = ImgInfo2.Model.PlaneParaInfo.SupOri;
tempSup( Img2IndexField( Img2IndexFieldOutliers)) = 0;
tempSup(:,end) = ImgInfo2.Model.PlaneParaInfo.SupOri(:,end);
tempSup(end,:) = ImgInfo2.Model.PlaneParaInfo.SupOri(end,:);
ImgInfo2.Model.PlaneParaInfo.SupOri = tempSup;

tempSup = ImgInfo2.Model.PlaneParaInfo.SupOri;
tempSup( Img2IndexTarget( Img2IndexTargetOutliers)) = 0;
tempSup(:,end) = ImgInfo2.Model.PlaneParaInfo.SupOri(:,end);
tempSup(end,:) = ImgInfo2.Model.PlaneParaInfo.SupOri(end,:);
ImgInfo2.Model.PlaneParaInfo.SupOri = tempSup;

tempSup = ImgInfo1.Model.PlaneParaInfo.SupOri;
tempSup( Img1IndexField( Img1IndexFieldOutliers)) = 0;
tempSup(:,end) = ImgInfo1.Model.PlaneParaInfo.SupOri(:,end);
tempSup(end,:) = ImgInfo1.Model.PlaneParaInfo.SupOri(end,:);
ImgInfo1.Model.PlaneParaInfo.SupOri = tempSup;


% assign model1 and model2
model1 = ImgInfo1.Model;
model2 = ImgInfo2.Model;

return;

