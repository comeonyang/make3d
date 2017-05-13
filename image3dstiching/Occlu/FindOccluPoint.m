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
function [POriReprojM TargetID FieldID FieldOccluPix OcluedDist FaceSetIDRemained] = FindOccluPoint(defaultPara, ScaleImg, ITarget, IField, F, ...
		TargetPointPix, FieldPointPix, INDTarget, INDField, ModelInfoTarget, ModelInfoField, Pair)

% Function find the occlu point 
% And output a region on epipolar line 

% initial parameters
FlagDisp = 1;
EpipolarThre = 0.1;
MaxFieldPointDist2Ray = 1;
[dump DepthScale(1) DepthScale(2)] = size(ModelInfoField.Model.PlaneParaInfo.Position3DFited);

NumTargetPointPix = size(TargetPointPix,2);
AllFieldPoint3D = ModelInfoField.Model.PlaneParaInfo.Position3DFited;
AllFieldPoint3D = AllFieldPoint3D(:,:);

% change FieldPoint3D into Target Coordinate
R_f_t = Pair.R';
T_f_t = -Pair.R'*Pair.T;
AllFieldPoint3D = R_f_t*AllFieldPoint3D + repmat( T_f_t, 1, size(AllFieldPoint3D,2)); % now in Target coordinate

% Initialize Variables
POriReprojM = [];
TargetID = [];
FieldID = [];
FaceSetIDRemained = [];
FieldOccluPix = [];
OcluedDist = [];

% Test Occlusion for each TargetPointPix point
fprintf(' Finding Occlusion Points ......');
FindOccluPointTime = tic;
for IndexTargetPointPix = 1:NumTargetPointPix
    
%     if IndexTargetPointPix/100 == round(IndexTargetPointPix/100)
%            IndexTargetPointPix
%     end
	% Initalize temp variables
	INDFieldCandidate = [];

    % 1) Pre-filtering the 2-d points in field image that are far from the epipolar line corresponding to the TargetPoint

	% Computing epipolar lines
	EpipolarLinePara = F*[ TargetPointPix(:,IndexTargetPointPix); 1];
	NormalizedEpipolarLinePara = EpipolarLinePara / norm(EpipolarLinePara(1:2));	

	% Find the RePojection of the Target Point in Field-Camera
	TargetPoint3D = ModelInfoTarget.Model.PlaneParaInfo.Position3DFited(:,INDTarget(:,IndexTargetPointPix));    % Size 3x1
	PointDepth(IndexTargetPointPix) = ModelInfoTarget.Model.PlaneParaInfo.FitDepth(INDTarget(:,IndexTargetPointPix));
        POriReproj =  defaultPara.InrinsicK2*(Pair.R*( TargetPoint3D) + Pair.T); 
        POriReproj = POriReproj(1:2,:)./POriReproj(3,:);

	% Modify the POriReproj to the closest point in the image
        % Since ReProjection of the Target Point might not be in size the Field image range
	POriReproj = Point2ImageRange(POriReproj, NormalizedEpipolarLinePara, ScaleImg); % some error/ Min fixed later//Solved 6/28

	% Finding Field points with distances to the epipolar line smaller than EpipolarThre*max(ScaleImg)
	Mask = abs(NormalizedEpipolarLinePara' *  [ FieldPointPix; ones(1, size(FieldPointPix,2))]) <= EpipolarThre*max(ScaleImg);
	INDFieldCandidate = INDField(Mask);
	FieldPoint3D = AllFieldPoint3D(:,INDFieldCandidate);	% 3 x NPoints % inFeild coordinate 
    
    
    % 2) Pre-filtering the 3-d points in field image that are far from the ray connecting target point to target image-origin.
   
 	% Find the distance of all the points from the target-origi ray
    	if ~isempty(FieldPoint3D)
        	u = TargetPoint3D' * FieldPoint3D / (TargetPoint3D' * TargetPoint3D);
        	distances = sqrt( sum( ( repmat(u,[3,1]) .* repmat(TargetPoint3D,[1,size(FieldPoint3D,2)]) - FieldPoint3D).^2 ,1));
        	occludingPointsMask = distances < MaxFieldPointDist2Ray;
		distances = distances(occludingPointsMask);
		[distances SortIND ] = sort(distances); 
		% update the new Field Candidate Points
		FieldPoint3D = FieldPoint3D(:,occludingPointsMask);
		FieldPoint3D = FieldPoint3D(:,SortIND);
		INDFieldCandidate = INDFieldCandidate(occludingPointsMask);
		INDFieldCandidate = INDFieldCandidate(SortIND);
    	end    
    
    % ============ Min's New trial =================
    % 3) Check triangals using FieldPoint3D are occluing the Ray or Not
	[ SubV SubH] = ind2sub(DepthScale, INDFieldCandidate);
%     length(INDFieldCandidate)
    i = 1;
    breakFlag = 0;
	while i <=length(INDFieldCandidate)	&& breakFlag == 0
%         i
		Ray =  TargetPoint3D./norm(TargetPoint3D);
		if (~(SubV(i) == 1 || SubH(i) == 1)	&& breakFlag ==0)
		% check top left square
			if ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i)-1, SubH(i)-1) == ...
				ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i), SubH(i))
				% upper right
				FaceSetID1 = INDFieldCandidate(i) - DepthScale(1) - 1;
				FaceSetID2 = INDFieldCandidate(i) - 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)-1 SubV(i)-1],[SubH(i) SubH(i)-1 SubH(i)]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
				% lower left
				FaceSetID1 = INDFieldCandidate(i) - DepthScale(1) - 1;
				FaceSetID2 = INDFieldCandidate(i) - DepthScale(1);
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)-1 SubV(i)],[SubH(i) SubH(i)-1 SubH(i)-1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			else
				% lower right
				FaceSetID1 = INDFieldCandidate(i) - DepthScale(1);
				FaceSetID2 = INDFieldCandidate(i) - 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i) SubV(i)-1],[SubH(i) SubH(i)-1 SubH(i)]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			end
		end

		if (~(SubV(i) == 1 || SubH(i) == DepthScale(2)) && breakFlag ==0)
		% check top right square
			if ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i)-1, SubH(i)) == ...
				ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i), SubH(i)+1)
				% lower left
				FaceSetID1 = INDFieldCandidate(i) - 1;
				FaceSetID2 = INDFieldCandidate(i) + DepthScale(1);
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)-1 SubV(i)],[SubH(i) SubH(i) SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			else
				% upper left
				FaceSetID1 = INDFieldCandidate(i) - 1;
				FaceSetID2 = INDFieldCandidate(i) + DepthScale(1) - 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)-1 SubV(i)-1],[SubH(i) SubH(i) SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
				% lower right
				FaceSetID1 = INDFieldCandidate(i) + DepthScale(1) - 1;
				FaceSetID2 = INDFieldCandidate(i) + DepthScale(1);
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)-1 SubV(i)],[SubH(i) SubH(i)+1 SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script

			end
		end
		
		if (~(SubV(i) == DepthScale(1) || SubH(i) == 1)	&& breakFlag ==0)
		% check bottom left square
			if ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i), SubH(i)-1) == ...
				ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i)+1, SubH(i))
				% upper right
				FaceSetID1 = INDFieldCandidate(i) - DepthScale(1);
				FaceSetID2 = INDFieldCandidate(i) + 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i) SubV(i)+1],[SubH(i) SubH(i)-1 SubH(i)]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script

			else
				% upper left
				FaceSetID1 = INDFieldCandidate(i) - DepthScale(1);
				FaceSetID2 = INDFieldCandidate(i) - DepthScale(1) + 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i) SubV(i)+1],[SubH(i) SubH(i)-1 SubH(i)-1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
				% lower right
				FaceSetID2 = INDFieldCandidate(i) - DepthScale(1) + 1;
				FaceSetID1 = INDFieldCandidate(i) + 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)+1 SubV(i)+1],[SubH(i) SubH(i)-1 SubH(i)]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			end
		end

		if (~(SubV(i) == DepthScale(1) || SubH(i) == DepthScale(2)) && breakFlag ==0)
		% check bottom left square
			if ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i), SubH(i)) == ...
				ModelInfoField.Model.PlaneParaInfo.SupOri(SubV(i)+1, SubH(i)+1)
				% upper right
				FaceSetID2 = INDFieldCandidate(i) + DepthScale(1) + 1;
				FaceSetID1 = INDFieldCandidate(i) + DepthScale(1);
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)+1 SubV(i)],[SubH(i) SubH(i)+1 SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
				% lower left
				FaceSetID2 = INDFieldCandidate(i) + 1;
				FaceSetID1 = INDFieldCandidate(i) + DepthScale(1) + 1;
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)+1 SubV(i)+1],[SubH(i) SubH(i) SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			else
				% upper left
				FaceSetID2 = INDFieldCandidate(i) + 1;
				FaceSetID1 = INDFieldCandidate(i) + DepthScale(1);
				TriInd = sub2ind(DepthScale,[ SubV(i) SubV(i)+1 SubV(i)],[SubH(i) SubH(i) SubH(i)+1]);
				[ CombinePara OccluPointDepth] = SolveOccluPoint( AllFieldPoint3D(:,TriInd),...
					Ray);
				OccluStorageScript; % min added script
			end
        end
%     i
    i = i + 1;
    end % end of second for loop
	% =====================================	
% IndexTargetPointPix
end % end of first for loop
disp([ '	' num2str( toc(FindOccluPointTime)) ' seconds']);
return;
