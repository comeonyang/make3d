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
function []= main( Fdir, PairList, Wrlname, defaultPara);

% This function is the main function tha handle:
% 1) Extracting info mation : Exif GPS IMU
% 2) Pair Image Metric Reconstruction

MainTime = tic; % Start counting time =============================
pause off;

if nargin < 2
    PairList = [];
end    
% initialize parameters
defaultPara.Wrlname = [defaultPara.Wrlname '_' Wrlname];

% initialize variables ============================================
ImgInfo = [];
[status currdir] = system([ 'ls ' Fdir '/data/ModelStatus/']);
if status
	[status currdir] = system([ 'mkdir ' Fdir '/data/ModelStatus/']); % Min add July 1st;
end
[status currdir] = system([ 'ls ' Fdir '/info/']);
if status
	[status currdir] = system([ 'mkdir ' Fdir '/info/']); % Min add Aug 17th;
end

if defaultPara.Flag.NewModelList
	ImgAddedList = cell(0); % cell array for names
else	
	load([ Fdir '/data/ModelStatus/' defaultPara.Wrlname '.mat']);
	NumImg = length(ModelStatus);
	for i = 1:NumImg
		ImgAddedList{i} = ModelStatus(i).ImgName;
	end
end

% initialize preacquired infomation ==============================
fprintf('Pre - Acquired Info .....');
if defaultPara.Flag.FlagImgInfoLoadPreStorage == 0
	[ImgInfo] = ExifExtractDir(Fdir,ImgInfo); % Exif info
    	% 1) Extracting info    
    	[ImgInfo] = ExtractGPSInfo(defaultPara, Fdir, ImgInfo); % GPS info
    	[ImgInfo] = ExtractRotationInfo( defaultPara, Fdir, ImgInfo); % IMU info
	% write into .kml file for google earth to display

	save([Fdir '/data/ImgInfo.mat'],'ImgInfo');
elseif defaultPara.Flag.FlagImgInfoLoadPreStorage == 1
	[ImgInfo] = ExifExtractDir(Fdir,ImgInfo); % Exif info
	% read from google earth .kml file
	temp = dir([Fdir '/info/*.kml']);
	[ ImgInfo ] = ExtractGoogleEarthInfo(ImgInfo, [Fdir '/info/' temp.name]);
    	[ ImgInfo ] = AllXWorld(defaultPara, ImgInfo); %X_world info
    	% [ImgInfo] = InfoInitialize(ImgInfo); % initialize everything
	save([Fdir '/data/ImgInfo.mat'],'ImgInfo');
elseif defaultPara.Flag.FlagImgInfoLoadPreStorage == 2
	[ImgInfo] = ExifExtractDir(Fdir,ImgInfo); % Exif info
    	[ImgInfo] = ExtractRotationInfo( defaultPara, Fdir, ImgInfo); % IMU info
	% read from google earth .kml file
	temp = dir([Fdir '/info/*.kml']);
	[ ImgInfo ] = ExtractGoogleEarthInfo(ImgInfo, [Fdir '/info/' temp.name]);
    	[ ImgInfo ] = AllXWorld(defaultPara, ImgInfo); %X_world info
    	% [ImgInfo] = InfoInitialize(ImgInfo); % initialize everything
	save([Fdir '/data/ImgInfo.mat'],'ImgInfo');
elseif defaultPara.Flag.FlagImgInfoLoadPreStorage == 3
	[ImgInfo] = ExifExtractDir(Fdir,ImgInfo); % Exif info
	save([Fdir '/data/ImgInfo.mat'],'ImgInfo');
else
    load([ Fdir '/data/ImgInfo.mat']); % // used for hand tune GPS and IMU info //Bad
end   
disp(['		' num2str( toc( MainTime)) ' seconds.']);

% SurfFeatures generation =========================================
fprintf('Calculating Surf Features ......');
MakeSureSurfDone(defaultPara, ImgInfo);
disp(['		' num2str( toc( MainTime)) ' seconds.']);

% start loop by getting input =====================================
NuPair = size(PairList,1);
ReadImgCount = 1;
while length(ImgAddedList) <= length(ImgInfo)
	
    % request user input
    if ReadImgCount > NuPair && defaultPara.Flag.NewInput
        disp('enter pair of image''s names');    
        if isempty(ImgA) || isempty(ImgB)    
            ImgA = input('add image A', 's');
            ImgB = input('add image B', 's');
        end
    elseif ReadImgCount <= NuPair
        ImgA = PairList{ReadImgCount,1};
        ImgB = PairList{ReadImgCount,2};
    else
	% if ReadImgCount > NuPair && ~defaultPara.Flag.NewInput 
	% ( Do not allow to add images manually)
        break;
    end
    
    if ReadImgCount == NuPair
        defaultPara.LastImgFlag = 1;
    end    
    ReadImgCount = ReadImgCount +1;
    
    % find ImgInfo index for ImgA ImgB
    [ImgAIndex] = ImgInfoIndexFromName(ImgInfo, ImgA);
    [ImgBIndex] = ImgInfoIndexFromName(ImgInfo, ImgB);

    % Metric reconstruction
    fprintf('Start Metric Reconstruction ....')
    [defaultPara NewImgInfo fail] = MetricRecon(defaultPara, ImgInfo([ImgAIndex ImgBIndex]));
    disp(['	' num2str( toc( MainTime)) ' seconds.']);
    if fail > 0
	disp('Failed ............... End of main.m');
	return;
    end
    ImgInfo(ImgAIndex ).Model = [];
    ImgInfo(ImgBIndex ).Model = [];   
    ImgInfo(ImgAIndex ) = NewImgInfo(1);
    ImgInfo(ImgBIndex ) = NewImgInfo(2);
	
    % add ImgAddedList with ImgA and ImgB
    ImgAddedList{end+1} = ImgA;
    ImgAddedList{end+1} = ImgB;
    
    % reset
    ImgA = [];
    ImgB = [];
    defaultPara.Flag.FlagFirstPair = 0;
end

pause on;
disp('End of main.m');
return;
