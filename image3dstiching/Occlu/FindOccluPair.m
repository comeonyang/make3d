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
function [PointPix1 PointDepth1 TargetIND1 POriReprojM1 FieldOccluPix1 OcluedDist1 OccluedFaceSetIND1 FieldFaceSetIDRemained1...
	  PointPix2 PointDepth2 TargetIND2 POriReprojM2 FieldOccluPix2 OcluedDist2 OccluedFaceSetIND2 FieldFaceSetIDRemained2]=...
                FindOccluPair(defaultPara, ModelInfo1, ModelInfo2, Pair, GlobalScale, SurfFlag)
            
% Function find the closed occlu point and output a region on epipolar line 
% to search for dense match
% Hacking Fix Parameters ==============================
Default.fy = 2400.2091651084;
Default.fx = 2407.3312729885838;
Default.Ox = 1110.7122391785729;%2272/2; %
Default.Oy = 833.72104535435108;%1704/2; %
Default.a_default = 2272/Default.fx; %0.70783777; %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
Default.b_default = 1704/Default.fy; %0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
Default.Ox_default = 1-Default.Ox/2272;%0.489272914; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;%0.488886982; % camera origin offset from the image center in vertical direction
[Default.VertYNuDepth Default.HoriXNuDepth] = size(ModelInfo1.Model.PlaneParaInfo.FitDepth);
% =====================================================

% loading Image
Img1 = ModelInfo1.ExifInfo.IDName;
I1 = imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']);
Img2 = ModelInfo2.ExifInfo.IDName;
I2 = imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']);
[Iy1 Ix1] = size(I1);
[DepthY1 DepthX1] = size(ModelInfo1.Model.Depth.FitDepth);
[Iy2 Ix2] = size(I2);
[DepthY2 DepthX2] = size(ModelInfo2.Model.Depth.FitDepth);

% Proper Scaling the Depth 3-d Position and PlaneParameters (all in global scale)
ModelInfo1.Model.PlaneParaInfo.PlanePara = ModelInfo1.Model.PlaneParaInfo.PlanePara/GlobalScale(1);
ModelInfo1.Model.PlaneParaInfo.FitDepth = ModelInfo1.Model.PlaneParaInfo.FitDepth*GlobalScale(1);
ModelInfo1.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo1.Model.PlaneParaInfo.Position3DFited*GlobalScale(1),[ 3 1 2]);
ModelInfo2.Model.PlaneParaInfo.PlanePara = ModelInfo2.Model.PlaneParaInfo.PlanePara/GlobalScale(2);
ModelInfo2.Model.PlaneParaInfo.FitDepth = ModelInfo2.Model.PlaneParaInfo.FitDepth*GlobalScale(1);
ModelInfo2.Model.PlaneParaInfo.Position3DFited = permute( ModelInfo2.Model.PlaneParaInfo.Position3DFited*GlobalScale(2),[ 3 1 2]);

% correct negative z positiom
ModelInfo1.Model.PlaneParaInfo.Position3DFited(3,:) = -ModelInfo1.Model.PlaneParaInfo.Position3DFited(3,:);
ModelInfo2.Model.PlaneParaInfo.Position3DFited(3,:) = -ModelInfo2.Model.PlaneParaInfo.Position3DFited(3,:);

% Define possible mathcing point pool
%  	SurfPoint (surfFeatures that haven't been matches)
%	FaceSetPoint (Point that compose the FaceSet in Wrl)
% Specify both Point in image Pixel position and their Index
% (Match index for Surf Features, DepthMap index for FaceSet)
    
% initalize FaceSetPoint---------------------
% ================Min to do: change to FaceSet base not Point base
FaceSetPickedIND1 = setdiff( 1:prod([DepthY1 DepthX1]), ...
			find(ModelInfo1.Model.PlaneParaInfo.SupOri == 0)); % used the ful index except the FaceSet Points that will not rendered
% ================================================================
Ray = ModelInfo1.Model.PlaneParaInfo.Ray;
RR =permute(Ray,[2 3 1]);
temp = RR(:,:,1:2)./repmat(RR(:,:,3),[1 1 2]);
FaceSetPoint1 = permute(temp./repmat(cat(3,Default.a_default,Default.b_default),[Default.VertYNuDepth Default.HoriXNuDepth 1])+...
		repmat(cat(3,Default.Ox_default,Default.Oy_default),[Default.VertYNuDepth Default.HoriXNuDepth 1]),[3 1 2]);
FaceSetPoint1(2,:) = 1- FaceSetPoint1(2,:);    
FaceSetPoint1 = FaceSetPoint1(:,:).*repmat([Ix1; Iy1],1,length(FaceSetPoint1(:,:)));
%FaceSetPoint1 = FaceSetPoint1(:,FaceSetPickedIND1);

% ================Min to do: change to FaceSet base not Point base (add July 4th)
FaceSetPickedIND2 = setdiff( 1:prod([DepthY2 DepthX2]), ...
			find(ModelInfo2.Model.PlaneParaInfo.SupOri == 0)); % used the ful index except the FaceSet Points that will not rendered
% ================================================================
Ray = ModelInfo2.Model.PlaneParaInfo.Ray;
RR =permute(Ray,[2 3 1]);
temp = RR(:,:,1:2)./repmat(RR(:,:,3),[1 1 2]);
FaceSetPoint2 = permute(temp./repmat(cat(3,Default.a_default,Default.b_default),[Default.VertYNuDepth Default.HoriXNuDepth 1])+...
		repmat(cat(3,Default.Ox_default,Default.Oy_default),[Default.VertYNuDepth Default.HoriXNuDepth 1]),[3 1 2]);
FaceSetPoint2(2,:) = 1- FaceSetPoint2(2,:);
FaceSetPoint2 = FaceSetPoint2(:,:).*repmat([Ix2; Iy2],1,length(FaceSetPoint2(:,:)));
%FaceSetPoint2 = FaceSetPoint2(:,FaceSetPickedIND2);
% ------------------------------------------------

% initialize SurfPoint---------------
if SurfFlag
	[SurfPoint1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % original features
	[SurfPoint2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % original features

	% Modified the ModelInfo to have the structure with the same number of SurfPoint
	ModelInfo1_surf = ModelInfo1;
	[SurfPoint2DepthMapApproxIND1] = ProjPosi2Mask( [Iy1 Ix1], [DepthY1 DepthX1], SurfPoint1);
	x_calib_1 = inv(defaultPara.InrinsicK1)*[ SurfPoint1; ones(1,size(SurfPoint1,2))];
	ray_1 = x_calib_1./repmat(sqrt( sum(x_calib_1.^2,1) ),3,1);
	Depth_1 = 1./sum( ModelInfo1.Model.PlaneParaInfo.PlanePara(:,ModelInfo1.Model.PlaneParaInfo.Sup2Para(ModelInfo1.Model.PlaneParaInfo.SupEpand(SurfPoint2DepthMapApproxIND1))).*ray_1, 1);
	Position3d_1 = ray_1.*repmat(Depth_1,3,1);
	% modified the structure of FitDepth and Position3DFited
	ModelInfo1_surf.Model.PlaneParaInfo.FitDepth = Depth_1;
	ModelInfo1_surf.Model.PlaneParaInfo.Position3DFited = Position3d_1;	

	ModelInfo2_surf = ModelInfo2;
    	[SurfPoint2DepthMapApproxIND2] = ProjPosi2Mask( [Iy2 Ix2], [DepthY2 DepthX2], SurfPoint2);
	x_calib_2 = inv(defaultPara.InrinsicK2)*[SurfPoint2; ones(1,length(SurfPoint2))];
	ray_2 = x_calib_2./repmat(sqrt( sum(x_calib_2.^2,1) ),3,1);
	Depth_2 = 1./sum( ModelInfo2.Model.PlaneParaInfo.PlanePara(:,ModelInfo2.Model.PlaneParaInfo.Sup2Para(ModelInfo2.Model.PlaneParaInfo.SupEpand(SurfPoint2DepthMapApproxIND2))).*ray_2, 1);
	Position3d_2 = ray_2.*repmat(Depth_2,3,1);
	% modified the structure of FitDepth and Position3DFited
	ModelInfo2_surf.Model.PlaneParaInfo.FitDepth = Depth_2;
	ModelInfo2_surf.Model.PlaneParaInfo.Position3DFited = Position3d_2;

	% get rid of matches SurfFeaturePoints (to aviod conflict in triangulation)
	SurfPickedIND1 = setdiff(1:size(SurfPoint1,2), Pair.matches(1,:)); % might use random pick surf Points if there are too many
	SurfPickedIND2 = setdiff(1:size(SurfPoint2,2), Pair.matches(2,:));
	%SurfPoint1 = SurfPoint1(:,SurfPickedIND1);
	%SurfPoint2 = SurfPoint2(:,SurfPickedIND2);
end

% Region by picking Img1
T1_hat = [[0 -Pair.T(3) Pair.T(2)];...
            [Pair.T(3) 0 -Pair.T(1)];...
            [-Pair.T(2) Pair.T(1) 0]];
F1 = inv(defaultPara.InrinsicK2)'*T1_hat*Pair.R*inv(defaultPara.InrinsicK1);
if SurfFlag
	[POriReprojM1 TargetIND1 OccluedFaceSetIND1 FieldOccluPix1 OcluedDist1 FieldFaceSetIDRemained1] = ...
		FindOccluPoint(defaultPara, [Iy2 Ix2], I1, I2, F1, ...
			SurfPoint1(:,SurfPickedIND1), FaceSetPoint2(:,FaceSetPickedIND2), SurfPickedIND1, FaceSetPickedIND2, ...
			ModelInfo1_surf, ModelInfo2, Pair);
	PointPix1 = SurfPoint1(:,TargetIND1);
	PointDepth1 = Depth_1(TargetIND1);
	% Structure Define:
	% POriReprojM1 TargetIND1 OccluedFaceSetIND1 FieldOccluPix1 OcluedDist1 
	% -- all of the size that single ray pass through both FaceSet 
	% OccluedFaceSetIND1 used in DepthMap size
	% "Notice" TargetIND1 is im surfFeature size
else
	[POriReprojM1 TargetIND1 OccluedFaceSetIND1 FieldOccluPix1 OcluedDist1 FieldFaceSetIDRemained1] = ...
		FindOccluPoint(defaultPara, [Iy2 Ix2], I1, I2, F1, ...
			FaceSetPoint1(:,FaceSetPickedIND1), FaceSetPoint2(:,FaceSetPickedIND2), FaceSetPickedIND1, FaceSetPickedIND2, ...
			ModelInfo1, ModelInfo2, Pair);
	PointPix1 = FaceSetPoint1(:,TargetIND1);
	PointDepth1 = ModelInfo1.Model.PlaneParaInfo.FitDepth(TargetIND1);
	% Structure Define:
	% POriReprojM1 TargetIND1 OccluedFaceSetIND1 FieldOccluPix1 OcluedDist1 
	% -- all of the size that single ray pass through both FaceSet 
	% OccluedFaceSetIND1 and TargetIND1 used in DepthMap size
end
save /afs/cs/group/reconstruction3d/scratch/testE/OccluPointFind.mat
% Region by picking Img2
Pair2.T = -Pair.R'*Pair.T;
Pair2.R = Pair.R';
Pair2.Xim = Pair.Xim([3 4 1 2],:);
if SurfFlag
	[POriReprojM2 TargetIND2 OccluedFaceSetIND2 FieldOccluPix2 OcluedDist2 FieldFaceSetIDRemained2] = ...
		FindOccluPoint(defaultPara, [Iy1 Ix1], I2, I1, F1', SurfPoint2(:,SurfPickedIND2), FaceSetPoint1(:,FaceSetPickedIND1), ...
			SurfPickedIND2, FaceSetPickedIND1, ModelInfo2_surf, ModelInfo1, Pair2);
	PointPix2 = SurfPoint2(:,TargetIND2);
	PointDepth2 = Depth_2(TargetIND2);
	% Structure Define:
	% POriReprojM2 TargetIND2 OccluedFaceSetIND2 FieldOccluPix2 OcluedDist2 
	% -- all of the size that single ray pass through both FaceSet 
	% OccluedFaceSetIND2 used in DepthMap size
	% "Notice" TargetIND2 is im surfFeature size
else
	[POriReprojM2 TargetIND2 OccluedFaceSetIND2 FieldOccluPix2 OcluedDist2 FieldFaceSetIDRemained2] = ...
		FindOccluPoint(defaultPara, [Iy1 Ix1], I2, I1, F1', FaceSetPoint2(:,FaceSetPickedIND2), FaceSetPoint1(:,FaceSetPickedIND1), ...
			FaceSetPickedIND2, FaceSetPickedIND1, ModelInfo2, ModelInfo1, Pair2);
	PointPix2 = FaceSetPoint2(:,TargetIND2);
	PointDepth2 = ModelInfo2.Model.PlaneParaInfo.FitDepth(TargetIND2);
	% Structure Define:
	%OccluedFaceSetIND2 used in DepthMap size
end
save /afs/cs/group/reconstruction3d/scratch/testE/OccluPointFind.mat


% Display
if false %defaultPara.Flag.FindOcclu
    MaskPositiveDist1 =  OccluDist1 >1;
    MaskPositiveDist2 =  OccluDist2 >1;  
    
    
    figure;
	MaxD1 = ones(1, sum(MaskPositiveDist1));
	MinD1 = MaxD1;
	MaxD2 = ones(1, sum(MaskPositiveDist2));
	MinD2 = MaxD2;
	
	dispMatchSearchRegin(I1, I2,[ PointPix1(:,MaskPositiveDist1); ones(1, sum(MaskPositiveDist1))], ...
        [ PointPix2(:,MaskPositiveDist2); ones(1, sum(MaskPositiveDist2))], ...
        Region1(:,MaskPositiveDist1), Region2(:,MaskPositiveDist2), F1, ...
        POriReprojM1(:,MaskPositiveDist1), MaxD1, PoccluM1(:,MaskPositiveDist1), MinD1, ...
        POriReprojM2(:,MaskPositiveDist2), MaxD2, PoccluM2(:,MaskPositiveDist2), MinD2, ...
        0, 'Stacking', 'v', 'Interactive', 0);
end
