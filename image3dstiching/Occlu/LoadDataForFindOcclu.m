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
function [ImgInfo1 ImgInfo2 Img1Index Img2Index Pair GlobalScale]=LoadDataForFindOcclu(defaultPara, WrlName, ImgName1, ImgName2, ...
						FlagPrestorage, PostFixStr)

% This function load every data FindOccluPair needs
if nargin < 6
	PostFixStr = 'NonMono';
end
% initalize the ImgInfo
ImgInfo = [];
if FlagPrestorage
	load([defaultPara.Fdir '/data/ImgInfo.mat']);
else
%	s = warning('off', 'all');
    [ImgInfo] = ExifExtractDir(defaultPara.Fdir,ImgInfo); % Exif info
    % 1) Extracting info
    [ImgInfo] = ExtractGPSInfo(defaultPara, defaultPara.Fdir, ImgInfo); % GPS info
    [ImgInfo] = ExtractRotationInfo( defaultPara, defaultPara.Fdir, ImgInfo); % IMU info
    [ ImgInfo ] = AllXWorld(defaultPara, ImgInfo); %X_world info
%	s = warning('on', 'all');
end
[Img1Index] = ImgInfoIndexFromName(ImgInfo, ImgName1);
[Img2Index] = ImgInfoIndexFromName(ImgInfo, ImgName2);
ImgInfo1 = ImgInfo(Img1Index);
ImgInfo2 = ImgInfo(Img2Index);
ImgInfo1.ExifInfo.IDName = ImgName1;
ImgInfo2.ExifInfo.IDName = ImgName2;

% load xxxx_NonMono.mat
load([ defaultPara.Fdir '/data/' ImgName1 '/' WrlName '_' ImgName1 '_' PostFixStr '.mat']);
ImgInfo1.Model = model;
load([ defaultPara.Fdir '/data/' ImgName2 '/' WrlName '_' ImgName2 '_' PostFixStr '.mat']);
ImgInfo2.Model = model;

% load xxx_xxx_PosrMatch.mat
% Notice the order
disp('Start loading Initial Pose info');
[Pair, Matches, Error] = LoadPoseMatch(defaultPara.Fdir, ImgName1, ImgName2);
Pair.matches = Matches;
if Error ~=0
   disp('Not a good Match data, might give bad R and T');
%    pause;
end

% load globalScale info in ModelStatus/xxx.mat
load([defaultPara.Fdir '/data/ModelStatus/' WrlName '.mat']);
GlobalScale =[];
Count = 0;
for i=1:length(ModelStatus)
	if strcmp( ModelStatus(i).ImgName, ImgName1)
		GlobalScale(1) = ModelStatus(i).Scale;
		R1 = ModelStatus(i).R;
		T1 = ModelStatus(i).T;
        Count = Count + 1;
	elseif strcmp( ModelStatus(i).ImgName, ImgName2)
		GlobalScale(2) = ModelStatus(i).Scale;
		R2 = ModelStatus(i).R;
                T2 = ModelStatus(i).T;
                Count = Count + 1;
    end
    
	if Count == 2
		break;
	end
end

% Min add Aug 22th =========
if ~exist('R2')	
    GlobalScale(2) = Pair.DepthScale(2)/Pair.DepthScale(1)*GlobalScale(1);
else
    % MainTain Global T nad R (coordinate 1 to 2)
    Pair.R = R2'*R1;
    Pair.T = (R2'*T1 - R2'*T2);
end    
% ==========================
return;
