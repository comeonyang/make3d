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
function [model1 model2] = RedundentPlanePick(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, ...
				Img1IndexTarget, Img2IndexField, Img2IndexTarget, Img1IndexField)

% [model1 model2] = RedundentPlanePick(defaultPara, ImgInfo1, ImgInfo2, ...
%                                Img1IndexTarget, Img2IndexField, Img2IndexTarget, Img1IndexField)

ErareTargetFlag = 1;
CheckOwnRayFlag = 1;
ErodeRemoveRegionFlag = 1;

% Setup the Position3D in Global Scale
Position3D1 = ImgInfo1.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo1.Model.PlaneParaInfo.Position3DFited*GlobalScale(1),[ 3 1 2]);
Position3D2 = ImgInfo2.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo1.Model.PlaneParaInfo.Position3DFited*GlobalScale(2),[ 3 1 2]);

% Pick out the Ray
Ray1 = ImgInfo1.Model.PlaneParaInfo.Ray(:,:);
Ray2 = ImgInfo2.Model.PlaneParaInfo.Ray(:,:);

% Calculate the Target Surface Norm
TargetFaceSetPoints1 = Position3D1(:,Img1IndexTarget(:));
TargerFaceSetPoints1 = reshape(TargerFaceSetPoints1, 3, 3, []);
TargerFaceSetTwoEdge1 = TargerFaceSetPoints1(:,1:2,:) - repmat( TargerFaceSetPoints1(:,3,:),[ 1 2 1]);
TargerSrufaceNormalVector1 = cross( TargerFaceSetTwoEdge1(:,1,:), TargerFaceSetTwoEdge1(:,2,:), 1);
TargerSrufaceNormalVector1 = permute( TargerSrufaceNormalVector1./repmat( sqrt( sum( TargerSrufaceNormalVector1.^2, 1)), [3 1 1]), [1 3 2]);% size 3 by NumberOfTarget

TargetFaceSetPoints2 = Position3D2(:,Img2IndexTarget(:));
TargerFaceSetPoints2 = reshape(TargerFaceSetPoints2, 3, 3, []);
TargerFaceSetTwoEdge2 = TargerFaceSetPoints2(:,1:2,:) - repmat( TargerFaceSetPoints2(:,3,:),[ 1 2 1]);
TargerSrufaceNormalVector2 = cross( TargerFaceSetTwoEdge2(:,1,:), TargerFaceSetTwoEdge2(:,2,:), 1);
TargerSrufaceNormalVector2 = permute( TargerSrufaceNormalVector2./repmat( sqrt( sum( TargerSrufaceNormalVector2.^2, 1)), [3 1 1]), [1 3 2]);% size 3 by NumberOfTarget

% Calculate the Field Surface Norm
FieldFaceSetPoints1 = Position3D2(:,Img2IndexField(:));
FieldFaceSetPoints1 = reshape(FieldFaceSetPoints1, 3, 3, []);
FieldFaceSetTwoEdge1 = FieldFaceSetPoints1(:,1:2,:) - repmat( FieldFaceSetPoints1(:,3,:),[ 1 2 1]);
FieldSrufaceNormalVector1 = cross( FieldFaceSetTwoEdge1(:,1,:), FieldFaceSetTwoEdge1(:,2,:), 1);
FieldSrufaceNormalVector1 = permute( FieldSrufaceNormalVector1./repmat( sqrt( sum( FieldSrufaceNormalVector1.^2, 1)), [3 1 1]), [1 3 2]);% size 3 by NumberOfTarget

FieldFaceSetPoints2 = Position3D1(:,Img1IndexField(:));
FieldFaceSetPoints2 = reshape(FieldFaceSetPoints2, 3, 3, []);
FieldFaceSetTwoEdge2 = FieldFaceSetPoints2(:,1:2,:) - repmat( FieldFaceSetPoints2(:,3,:),[ 1 2 1]);
FieldSrufaceNormalVector2 = cross( FieldFaceSetTwoEdge2(:,1,:), FieldFaceSetTwoEdge2(:,2,:), 1);
FieldSrufaceNormalVector2 = permute( FieldSrufaceNormalVector2./repmat( sqrt( sum( FieldSrufaceNormalVector2.^2, 1)), [3 1 1]), [1 3 2]);% size 3 by NumberOfTarget

% Calculate the portion that align with the ray
SurfNorm1Target2Ray = abs( sum(TargerSrufaceNormalVector1.*Ray1(:,Img1IndexTarget), 1) );
if CheckOwnRayFlag
	SurfNorm1Field2Ray = abs( sum(FieldSrufaceNormalVector1.*Ray2(:,Img2IndexField), 1) );
else
	SurfNorm1Field2Ray = abs( sum( ( Pair.R*FieldSrufaceNormalVector1 + repmat(Pair.T, 1, size(FieldSrufaceNormalVector1, 2) )).*Ray1(:,Img1IndexTarget), 1) );
end

% Calculate the portion that align with the ray
SurfNorm2Target2Ray =  sum( TargerSrufaceNormalVector2.*Ray2(:,Img2IndexTarget), 1) ;
if CheckOwnRayFlag
	SurfNorm2Field2Ray = abs( sum( FieldSrufaceNormalVector2.*Ray1(:,Img1IndexField), 1) );
else
	SurfNorm2Field2Ray = abs( sum( (Pair.R'*FieldSrufaceNormalVector2 + repmat(-Pair.R'*Pair.T, 1, size(FieldSrufaceNormalVector2, 2) )).*Ray2(:,Img2IndexTarget), 1) );
end

% ===============================================================================

% Pick out the Plane has a bigger normalized PlanPrameter of the Ray
if ErareTargetFlag
	Img1IndexTargetOutliers = abs(SurfNorm1Target2Ray) < abs(SurfNorm1Field2Ray); 
else
	Img1IndexTargetOutliers = Inf < abs(SurfNorm1Field2Ray);
end
Img2IndexFieldOutliers = ~Img1IndexTargetOutliers;

if ErareTargetFlag
	Img2IndexTargetOutliers = abs(SurfNorm2Target2Ray) < abs(SurfNorm2Field2Ray); 
else
	Img2IndexTargetOutliers = Inf < abs(SurfNorm2Field2Ray); 
end
Img1IndexFieldOutliers = ~Img2IndexTargetOutliers;


if ErodeRemoveRegionFlag
   OriSup = ones(size(ImgInfo1.Model.PlaneParaInfo.SupOri));
   OriSup( Img1IndexTarget( Img1IndexTargetOutliers)) = 0;
   se = strel('ball',5,5);
   OriSup = imdilate(OriSup,se);
   Img1IndexTargetNew = find(OriSup == 0);
   
   OriSup = ones(size(ImgInfo2.Model.PlaneParaInfo.SupOri));
   OriSup(  Img2IndexField( Img2IndexFieldOutliers)) = 0;
   se = strel('ball',5,5);
   OriSup = imdilate(OriSup,se);
   Img2IndexFieldNew = find(OriSup == 0);
   
   OriSup = ones(size(ImgInfo2.Model.PlaneParaInfo.SupOri));
   OriSup(  Img2IndexTarget( Img2IndexTargetOutliers)) = 0;
   se = strel('ball',5,5);
   OriSup = imdilate(OriSup,se);
   Img2IndexTargetNew = find(OriSup == 0);
   
   OriSup = ones(size(ImgInfo1.Model.PlaneParaInfo.SupOri));
   OriSup( Img1IndexField( Img1IndexFieldOutliers)) = 0;
   se = strel('ball',5,5);
   OriSup = imdilate(OriSup,se);
   Img1IndexFieldNew = find(OriSup == 0);
end    

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

