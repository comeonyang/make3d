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
function [ImgInfo] = ReInference(defaultPara, R, T, Xim, Model, TriDepth, GroundLevel, CoordinateFromRef, FlagRender)

% This function run the inference in two phase
%	Phase 1: Propagate triangulate points info
%	Phase 2: Refine the non-triangulated Sup to Vertical and Horizontal setting
%		 and the multiple image overlap information

% Input:
%	R - rotaion from coordinate ImgInfo(2) to ImgInfo(1)
%	T - translation from coordinate ImgInfo(2) to ImgInfo(1)
%	Xim - matches in pixel coordinate
%	Model - depth previous constrains
%	TriDepth - triangulated info
%	GroundLevel - ground level info (force ground to have the same level)
%	CoordinateFromRef - Rotaion and translation ( [3x3 3x1]) from ImgInfo(1) to the world reference (the first image)a
%	FlagRender - rendering the vrml or not
% Return:
%	ImgInfo - add constrain;
%	

% initialize parameters -----------------------------------------------------
scale = 1; 
StickHori = 5;  %0.1; % sticking power in horizontal direction
StickVert = 5;     % sticking power in vertical direction
Center = 10; % Co-Planar weight at the Center of each superpixel
TriangulatedWeight = 20;
GroundLevelWright = 500;
GroundWeight = 20;
VertWeight = 10
HoriConf = 1; % set the confidant of the learned depth at the middle in Horizontal direction of the image
VertConf = 0.01; % set the confidant of the learned depth at the top of the image
[ VertYNuDepth HoriXNuDepth] = size(Model(1).Depth.FitDepth);
mapVert = linspace(VertConf,1,VertYNuDepth); % modeling the gravity prior
mapHori = [linspace(HoriConf,1,round(HoriXNuDepth/2)) fliplr(linspace(HoriConf,1,HoriXNuDepth-round(HoriXNuDepth/2)))];
% ========set the range of depth that our model in 
ClosestDist = defaultPara.Closestdist;
FarestDist = defaultPara.FarestDist; % //Min: do not need it since depths are all been processed
% ================================================
ceiling = 0*VertYNuDepth; % set the position of the ceiling, related to No plane coming back constrain % changed for newchurch
if ~isempty(Model(1).MultiScaleSupTable)
   MultiScaleFlag = true;
   WeiV = 2*ones(1,size(Model(1).MultiScaleSupTable,2)-1);
else
   MultiScaleFlag = false;
   WeiV = 1;
end
WeiV(1,1:2:end) = 6; % emphasize the middle scale three times smaller than large scale
WeiV = WeiV./sum(WeiV);% normalize if pair of superpixels have same index in all the scale, their weight will be 10
ShiftStick = -.1;  % between -1 and 0, more means more smoothing.
ShiftCoP = -.5;  % between -1 and 0, more means more smoothing.
gravity =true; % if true, apply the HoriConf and VertConf linear scale weight
% =======================================
groundThreshold = cos([ zeros(1, VertYNuDepth - ceil(VertYNuDepth/2)+10) ...
                                linspace(0,15,ceil(VertYNuDepth/2)-10)]*pi/180);
    %  v1 15 v2 20 too big v3 20 to ensure non misclassified as ground.
    %  verticalThreshold = cos(linspace(5,55,Default.VertYNuDepth)*pi/180); % give a vector of size 55 in top to down : 
verticalThreshold = cos([ 5*ones(1,VertYNuDepth - ceil(VertYNuDepth/2)) ...
                                linspace(5,55,ceil(VertYNuDepth/2))]*pi/180);
        % give a vector of size 55 in top to down : 
        % 50 means suface norm away from y axis more than 50 degree
% ===========================================================================================================================================
if strcmp(defaultPara.InitialDepth,'FitDepth')
	CleanedDepthMap = Model(1).Depth.FitDepth;
else
	CleanedDepthMap = Model(1).Depth.RawDepth;
end

% -----------------------------------------------------------------------------------

% Pre-Processing --------------------------------------------------------------------
% Clean the Sup near sky
maskSky = Model(1).Sup == 0;
maskSkyEroded = imerode(maskSky, strel('disk', 4) );
SupEpand = ExpandSup2Sky(Sup,maskSkyEroded);
NuPatch = HoriXNuDepth*VertYNuDepth-sum(maskSky(:));
NuSup = setdiff(unique(Sup)',0);
NuSup = sort(NuSup);
NuSupSize = size(NuSup,2);
% Sup index and planeParameter index inverse map
Sup2Para = sparse(1,max(Sup(:)));
Sup2Para(NuSup) = 1:NuSupSize;

Posi3D = im_cr2w_cr(Model(1).Depth.FitDepth, permute(Model(1).Ray,[2 3 1])); % (3 by VertNuDepth HoriNuDepth)
% -----------------------------------------------------------------------------------

% Generate the Matrix for MRF -------------------------------------------------------
CleanedDepthMap = Model(1).Depth.FitDepth; 
PosiM = sparse(0,0); % Position matrix: first self term objective
RayAllM = sparse(0,0); % all the Ray for regular grid to set depth constrain
ScalingTerm = sparse( 0, 1);
YPointer = [];
YPosition = [];
beta = [];
for i = NuSup
    mask = SupEpand ==i; % include the Ray that will be use to expand the NonSky
    RayAllM = blkdiag( RayAllM, Ray(:,mask)');
    mask = Sup ==i; % Not include the Ray that will be use to expand the NonSky    
    [yt x] = find(mask);
    CenterX = round(median(x));
    CenterY = round(median(yt));
    YPointer = [YPointer; CenterY >= ceiling]; % Y is zero in the top ceiling default to be 0 as the top row in the image
    YPosition = [YPosition; CenterY];
    mask(isnan(CleanedDepthMap)) = false;
    SupNuPatch(i) = sum(mask(:));
    % find center point
    [yt x] = find(mask);
    CenterX = round(median(x));
    CenterY = round(median(yt));

  if ~all(mask(:)==0)
    if gravity
      if any(CleanedDepthMap(mask) <=0)
         CleanedDepthMap(mask)
      end
      PosiM = blkdiag(PosiM,Posi3D(:,mask)');%.*repmat( mapVert(yt)',[1 3]).*repmat( mapHori(x)',[1 3]));
      if SupMatched == i
         ScalingTerm = [ScalingTerm; ones( sum(mask(:)), 1)];
      else
         ScalingTerm = [ScalingTerm; zeros( sum(mask(:)), 1)];
      end
    else
      PosiM = blkdiag(PosiM,Posi3D(:,mask)');
      if SupMatched == i
         ScalingTerm = [ScalingTerm; ones( sum(mask(:)), 1)];
      else
         ScalingTerm = [ScalingTerm; zeros( sum(mask(:)), 1)];
      end
    end
  else
     PosiM = blkdiag(PosiM, Posi3D(:,mask)');
  end
end
YPointer(YPointer==0) = -1;

% =================Building up the triangulated  and sampled ground constrain ==========================
PosiTriangulatedM = sparse( 0, 3*NuSupSize);
count = 1;
for i = SupMatched'
    temp = sparse(1, 3*NuSupSize);
    if Sup2Para(i)*3>3*NuSupSize || Sup2Para(i)*3-2<0
       Sup2Para(i)
       NuSupSize
       count = count + 1;
       continue;
    end
%     Sup2Para(i)
%     i
    temp( (Sup2Para(i)*3-2):(Sup2Para(i)*3) ) = RayMatched(count,:)*ClosestDepth(count);
    PosiTriangulatedM = [PosiTriangulatedM; temp];
    count = count + 1;
end
PosiSampledGroundM = sparse( 0, 3*NuSupSize);
count = 1;
for i = SampledGroundSupMatched'
    temp = sparse(1, 3*NuSupSize);
    if (Sup2Para(i)*3)>(3*NuSupSize) || (Sup2Para(i)*3-2)<0
       Sup2Para(i)
       NuSupSize
    end
    temp( (Sup2Para(i)*3-2):(Sup2Para(i)*3) ) = SampledGroundRayMatched(count,:)*SampledGroundClosestDepth(count);
    PosiSampledGroundM = [PosiSampledGroundM; temp];
    count = count + 1;
end

% ================================================================================================


