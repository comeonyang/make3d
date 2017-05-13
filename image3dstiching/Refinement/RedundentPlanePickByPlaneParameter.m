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

% take out the Plane Parameter to check the orientationa==========================
PlanePara1Target = ImgInfo1.Model.PlaneParaInfo.PlanePara(:,ImgInfo1.Model.PlaneParaInfo.Sup2Para(...
                   ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexTarget)));
PlanePara1Field = ImgInfo2.Model.PlaneParaInfo.PlanePara(:,ImgInfo2.Model.PlaneParaInfo.Sup2Para(...
                   ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexField)));
% correct Scale to Global Scale
PlanePara1Target = PlanePara1Target/GlobalScale(1);
PlanePara1Field = PlanePara1Field/GlobalScale(2);
% Transfer Feild Plane Parameters into Target coordinate
PlanePara1Field_Ori = PlanePara1Field;
PlanePara1Field = Pair.R'*PlanePara1Field./repmat(1- Pair.T'*PlanePara1Field, 3, 1);

% normalized to unit vector
PlanePara1Target = PlanePara1Target./repmat(sqrt( sum( PlanePara1Target.^2,1)), 3, 1);
PlanePara1Field_Ori = PlanePara1Field_Ori./repmat(sqrt( sum( PlanePara1Field_Ori.^2,1)), 3, 1);
PlanePara1Field = PlanePara1Field./repmat(sqrt( sum( PlanePara1Field.^2,1)), 3, 1);

% Calculate the portion that align with the ray
Ray1 = ImgInfo1.Model.PlaneParaInfo.Ray(:,:);
PlanePara1Target2Ray = abs( sum(PlanePara1Target.*Ray1(:,Img1IndexTarget), 1) );
if CheckOwnRayFlag
	Ray1Field = ImgInfo2.Model.PlaneParaInfo.Ray(:,:);
	PlanePara1Field2Ray = abs( sum(PlanePara1Field_Ori.*Ray1Field(:,Img2IndexField), 1) );
else
	PlanePara1Field2Ray = abs( sum(PlanePara1Field.*Ray1(:,Img1IndexTarget), 1) );
end


PlanePara2Target = ImgInfo2.Model.PlaneParaInfo.PlanePara(:,ImgInfo2.Model.PlaneParaInfo.Sup2Para(...
                    ImgInfo2.Model.PlaneParaInfo.SupEpand(Img2IndexTarget)));
PlanePara2Field = ImgInfo1.Model.PlaneParaInfo.PlanePara(:,ImgInfo1.Model.PlaneParaInfo.Sup2Para(...
                    ImgInfo1.Model.PlaneParaInfo.SupEpand(Img1IndexField)));

% correct Scale to Global Scale
PlanePara2Target = PlanePara2Target/GlobalScale(2);
PlanePara2Field = PlanePara2Field/GlobalScale(1);
% Transfer Feild Plane Parameters into Target coordinate
PlanePara2Field_Ori = PlanePara2Field;
PlanePara2Field = Pair.R*PlanePara2Field./repmat(1- ( -Pair.R'*Pair.T )'*PlanePara2Field, 3, 1);

% normalized to unit vector
PlanePara2Target = PlanePara2Target./repmat(sqrt( sum( PlanePara2Target.^2,1)), 3, 1);
PlanePara2Field_Ori = PlanePara2Field_Ori./repmat(sqrt( sum( PlanePara2Field_Ori.^2,1)), 3, 1);
PlanePara2Field = PlanePara2Field./repmat(sqrt( sum( PlanePara2Field.^2,1)), 3, 1);

% Calculate the portion that align with the ray
Ray2 = ImgInfo2.Model.PlaneParaInfo.Ray(:,:);
PlanePara2Target2Ray =  sum(PlanePara2Target.*Ray2(:,Img2IndexTarget), 1) ;
if CheckOwnRayFlag
	Ray2Field = ImgInfo1.Model.PlaneParaInfo.Ray(:,:);
	PlanePara2Field2Ray = abs( sum(PlanePara2Field_Ori.*Ray2Field(:,Img1IndexField), 1) );
else
	PlanePara2Field2Ray =  sum(PlanePara2Field.*Ray2(:,Img2IndexTarget), 1) ;
end

% ===============================================================================

% Pick out the Plane has a bigger normalized PlanPrameter of the Ray
if ErareTargetFlag
	Img1IndexTargetOutliers = abs(PlanePara1Target2Ray) < abs(PlanePara1Field2Ray); 
else
	Img1IndexTargetOutliers = Inf < abs(PlanePara1Field2Ray);
end
Img2IndexFieldOutliers = ~Img1IndexTargetOutliers;

if ErareTargetFlag
	Img2IndexTargetOutliers = abs(PlanePara2Target2Ray) < abs(PlanePara2Field2Ray); 
else
	Img2IndexTargetOutliers = Inf < abs(PlanePara2Field2Ray); 
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

