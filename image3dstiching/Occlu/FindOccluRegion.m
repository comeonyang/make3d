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
function [Region POriReprojM PoccluM OccluDist] = FindOccluRegion(defaultPara, ScaleImg, ITarget, F, TargetPointPix, FieldPointPix, INDTarget, INDField, ModelInfoTarget, ModelInfoField, Pair)

% Function find the closed occlu point and output a region on epipolar line 
% to search for dense match

% find features close to epipoaeline

% initial parameters
FlagDisp = 1;
EpipolarThre = 0.1;
ScaleSlab = 1;
VertSearchSlab = 0.1;
HoriSearchSlab = 0.05;

% Initialize Variables
Region = [];
POriReprojM= [];
PoccluM = [];
OccluDist = [];

% Generate a region for each TargetPointPix point
NumTargetPointPix = size(TargetPointPix,2);
for IndexTargetPointPix = 1:NumTargetPointPix
% iterating over all points in target image

	% Pre-filtering the 3-d points in field image that are far from the ray connecting target point to target image-origin.
	%closeToRay = true;
	% depth
        %[DTargetPick INDTargetPick] = PorjPosi2Depth([Iy Ix], [DepthY DepthX], TargetPointPix(:,IndexTargetPointPix), 
	%								ModelInfoTarget.Model.PlaneParaInfo.FitDepth);
	% computing epipolar lines
	EpipolarLinePara = F*[ TargetPointPix(:,IndexTargetPointPix); 1];
	EpipolarLinePara = EpipolarLinePara / norm(EpipolarLinePara(1:2));	
	TargetPoint3D = ModelInfoTarget.Model.PlaneParaInfo.Position3DFited(:,INDTarget(:,IndexTargetPointPix));    % 3x1

	%FieldPoint3D = ModelInfoField.Model.PlaneParaInfo.Position3DFitted;	% 55x305x3
	
	% finding the distance to the line

	% Find the RePojection of the TargetPick Point in Field-Camera
        x_calib = inv(defaultPara.InrinsicK1)*[ TargetPointPix(:,IndexTargetPointPix); 1];
        ray = x_calib/norm(x_calib);
        POriReproj =  defaultPara.InrinsicK2*(Pair.R*( ray*ModelInfoTarget.Model.PlaneParaInfo.FitDepth( INDTarget(:,IndexTargetPointPix) )) + Pair.T); 
%         ModelInfoTarget.Model.PlaneParaInfo.Position3DFited(:,INDTarget(:
%         ,IndexTargetPointPix))
        POriReproj = POriReproj(1:2,:)./POriReproj(3,:);
	% modify the POriReproj to the closest point in the image
	% First check Horizontal then check Vertical
	LineVector = [ EpipolarLinePara(2) -EpipolarLinePara(1)];
	if POriReproj(2) > ScaleImg(2) % make POriReproj(2) = ScaleImg(2)
		shift = (ScaleImg(2) - POriReproj(2))/LineVector(2);
		POriReproj(1) = POriReproj(1) + shift*LineVector(1);
		POriReproj(2) = ScaleImg(2);
	elseif POriReproj(2) < 1
		shift = (1 - POriReproj(2))/LineVector(2);
		POriReproj(1) = POriReproj(1) + shift*LineVector(1);
		POriReproj(2) = 1;
	end
	
	if POriReproj(1) > ScaleImg(1) % make POriReproj(1) = ScaleImg(1)
                shift = (ScaleImg(1) - POriReproj(1))/LineVector(1);
                POriReproj(2) = POriReproj(2) + shift*LineVector(2);
                POriReproj(1) = ScaleImg(1);
        elseif POriReproj(1) < 1
                shift = (1 - POriReproj(1))/LineVector(1);
                POriReproj(2) = POriReproj(2) + shift*LineVector(2);
                POriReproj(1) = 1;
        end

	if abs(EpipolarLinePara' * [POriReproj ; ones(1, 1)]) > 0.1 
        	disp('error of POriReproj');
	end
		
	% Find the ReProjection of the origion of Target camer in Field camera
	%if Pair.T(3) < 0 

	%else

	%end

	Mask = abs(EpipolarLinePara' *  [ FieldPointPix; ones(1, size(FieldPointPix,2))]) <= EpipolarThre*max(ScaleImg);

    FieldPoint3D = ModelInfoField.Model.PlaneParaInfo.Position3DFited(:,INDField(Mask));	% 3 x NPoints % inFeild coordinate
	FieldPointPixMask = FieldPointPix(:, Mask);
	% change FieldPoint3D into Target Coordinate
	R_f_t = Pair.R';
	T_f_t = -Pair.R'*Pair.T;
	FieldPoint3D = R_f_t*FieldPoint3D + repmat( T_f_t, 1, size(FieldPoint3D,2)); % now in Target coordinate
    
    % find if the points that occlude
    
    % first find the distance of all the point from the target-origi ray
    if ~isempty(FieldPoint3D)
        u = TargetPoint3D' * FieldPoint3D / (TargetPoint3D' * TargetPoint3D);
        distances = sqrt( sum( ( repmat(u,[3,1]) .* repmat(TargetPoint3D,[1,size(FieldPoint3D,2)]) - FieldPoint3D).^2 ,1));
        [ V occludingPointsMask] = min(distances);
%         occludingPointsMask = distances < 0.1;    % MAGIC number. 1 meter in world units
    else
        occludingPointsMask = [];
    end    
    
    % ============
    % check the plan parameter to make sure alpha'*x - 1 > 0
    % check is the occlusion reprojection has the same planeparameter index
    % ================
    
    if ~isempty(occludingPointsMask) && sum(occludingPointsMask) ~=0% occludsion case
        candidateFieldPoints3D = FieldPointPixMask(:,occludingPointsMask);
        candidateU = u(occludingPointsMask); % fraction of TargetPoint3D
        [tmp, IND] = min(candidateU);
        tmp = tmp*norm(TargetPoint3D); % change to real unit meters
    %  if tmp < ModelInfoTarget.Model.PlaneParaInfo.FitDepth(INDTarget(:,IndexTargetPointPix))
    
	% find the closest occlusion to camera-Target origin.
	%[Pocclu] = ClosetOcclusion( defaultPara, ScaleImg, Pair, ModelInfoTarget, ...
	%			TargetPointPix(:,IndexTargetPointPix), INDTarget(:,IndexTargetPointPix), ModelInfoField, FieldPointPix(:, Mask), INDField(Mask));

	%if closeToRay
		% consider intersection with the field superpixel of that point
	%end
	
    %if ~isempty(Pocclu)
        % DTargetPick and Pair.T both are rescaled to global scale
	
        % Define rectangular region form Pocclu and POriReproj
        PoccluNearBy = candidateFieldPoints3D(:,IND);
        if isempty(PoccluNearBy)            
            disp('error')
        end
	% need to Porject PoccluNearBy to the closest position on the epipolar line
	Pocclu = - ( EpipolarLinePara' *  [ PoccluNearBy; ones(1, 1)])* EpipolarLinePara(1:2) + PoccluNearBy;
    	if Pocclu(2) > ScaleImg(2) % make Pocclu(2) = ScaleImg(2)
		shift = (ScaleImg(2) - Pocclu(2))/LineVector(2);
		Pocclu(1) = Pocclu(1) + shift*LineVector(1);
		Pocclu(2) = ScaleImg(2);
	elseif Pocclu(2) < 1
		shift = (1 - Pocclu(2))/LineVector(2);
		Pocclu(1) = Pocclu(1) + shift*LineVector(1);
		Pocclu(2) = 1;
	end
	
	if Pocclu(1) > ScaleImg(1) % make Pocclu(1) = ScaleImg(1)
                shift = (ScaleImg(1) - Pocclu(1))/LineVector(1);
                Pocclu(2) = Pocclu(2) + shift*LineVector(2);
                Pocclu(1) = ScaleImg(1);
	elseif Pocclu(1) < 1
                shift = (1 - Pocclu(1))/LineVector(1);
                Pocclu(2) = Pocclu(2) + shift*LineVector(2);
                Pocclu(1) = 1;
    	end
    
	if abs(EpipolarLinePara' * [Pocclu ; ones(1, 1)]) > 0.1
        	disp('error of Pocclu');
	end

        vari = [[ 0 0 HoriSearchSlab -HoriSearchSlab];...
        	[VertSearchSlab -VertSearchSlab 0 0]];	
        vari = vari*max(ScaleImg)*ScaleSlab;
        MaxV = min(max([Pocclu(2)+vari(2,:) POriReproj(2)+vari(2,:)]),ScaleImg(1));
        MinV = max(min([Pocclu(2)+vari(2,:) POriReproj(2)+vari(2,:)]),1);
        MaxH = min(max([Pocclu(1)+vari(1,:) POriReproj(1)+vari(1,:)]),ScaleImg(2));
        MinH = max(min([Pocclu(1)+vari(1,:) POriReproj(1)+vari(1,:)]),1);
        %Region = [Region [MaxH; MinH; MaxV; MinV]];
        Region = [Region [MinH; MaxH; MinV; MaxV]];
        PoccluM = [PoccluM Pocclu];
        POriReprojM = [POriReprojM POriReproj];
        OccluDist =[ OccluDist (ModelInfoTarget.Model.PlaneParaInfo.FitDepth(INDTarget(:,IndexTargetPointPix)) -tmp)];
     % else
     %   Region = [Region zeros(4,1)];
     %   PoccluM = [PoccluM zeros(2,1)];
     %   POriReprojM = [POriReprojM zeros(2,1)];
     %   OccluDist = [ OccluDist zeros(1,1)];
     % end  
    else % case of no occlusion
        Region = [Region zeros(4,1)];
        PoccluM = [PoccluM zeros(2,1)];
        POriReprojM = [POriReprojM zeros(2,1)];
        OccluDist = [ OccluDist zeros(1,1)];
    end
end

return;
