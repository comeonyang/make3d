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
function SingleImgReInference(defaultPara, model, Img, GroundLevel, Wrlname, R, T, Scale,PostFixStr)

% This funciton run the inference by loading all constrain in the model variable

% Fix parameters
% ==== General Parameters ===============
Default.RenderFlag =  defaultPara.RenderFlag;
[Default.VertYNuDepth Default.HoriXNuDepth] = size(model.Depth.FitDepth);
Default.fy = 2400.2091651084;
Default.fx = 2407.3312729885838;
Default.Ox = 1110.7122391785729;%2272/2; %
Default.Oy = 833.72104535435108;%1704/2; %
Default.a_default = 2272/Default.fx; %0.70783777; %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
Default.b_default = 1704/Default.fy; %0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
Default.Ox_default = 1-Default.Ox/2272;%0.489272914; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;%0.488886982; % camera origin offset from the image center in vertical direction
Default.MinTriEffectPercent = 10; % 5 % higher the CoPlaner term % change to 10 make sure tri area are flat not jaggy
Default.FarestTriDist = 1;%0.1; % unit in pixel

% TriCountSupThre depends on the total number of model.Constrain.RayMatche
Default.TriCountSupThre = size(model.Constrain.RayMatche,1)*0.1;
% ========================================== First Model

    % highly changable
    Default.OutPutFolder = [ defaultPara.OutPutFolder Img '/'];
    Default.ScratchFolder = Default.OutPutFolder;
    Default.Wrlname{1} = [Wrlname '_' Img '_' PostFixStr];
    Default.Flag.AfterInferenceStorage = 0;

    AappendOpt = 0; % not append since generate seperate .wrl files

	% simply constrain of the PairList matches
    Aconstrain.RayMatched = model.Constrain.RayMatche;
    Aconstrain.Depth_modified = model.Constrain.Depth_modified;
    Aconstrain.SupMatched = model.Constrain.SupMatched;

	% new constraian of the occlusion Matches
    Aconstrain.OccluRayMatched = cell2mat(model.ConstrainOccluMatch.RayMatche(:));
    Aconstrain.OccluDepth_modified = cell2mat(model.ConstrainOccluMatch.Depth_modified(:)');
    Aconstrain.OccluSupMatched = cell2mat(model.ConstrainOccluMatch.SupMatched(:));
	
    ASup = model.Sup;
    ASupOri = ASup;
    AdepthMap = model.Depth.RawDepth; % or RawDepth
    ARayOri = model.Ray;
    ARayAll = ARayOri;
    ASupNeighborTable = model.SupNeighborTable;
    AmaskSky = model.maskSky;
    AmaskG = model.maskG;
    AMultiScaleSupTable = model.MultiScaleSupTable;
    
    Default.filename{1} = ['../' Img];
% Calling PairReInferenceSepRender

[ model.PlaneParaInfo] = PlaneParaMRFTriangulateOneShot( Default, R, T, AappendOpt, ...
                           [ Aconstrain.RayMatched; Aconstrain.OccluRayMatched],...
                           [ Aconstrain.Depth_modified Aconstrain.OccluDepth_modified]',...
                           [ Aconstrain.SupMatched; Aconstrain.OccluSupMatched],...
                           [ ], [ ], [ ], [],...
                           ASup, ASupOri, [], AdepthMap, zeros(size(AdepthMap)), ARayOri, ARayAll, ...
                           ASupNeighborTable, [], AmaskSky, AmaskG,...
                           'cvx_allL1Norm',1,...
                           [], [], AMultiScaleSupTable, [], [], [], false, Scale, R, GroundLevel);% eye(3), GroundLevel might be wrong //Min check

    % Important Storage the Triangulated info in local scale
    save( [defaultPara.ScratchFolder Img '/' Wrlname '_' Img '_' PostFixStr '.mat'], 'model'); % add prefix defaultPara.Wrlname to distinguish model
return;
