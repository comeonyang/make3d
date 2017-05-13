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
function makeDepthImagesSph_ImgbyImg(dataNumber, displayFlag)
% global depthBins;
% global countBuf;
if nargin < 2
    displayFlag = false;
end
if nargin < 1
    dataNumber = 1;
end
%% we first read the log file and then read the depth values in a
%% two-dimensional array after binning and averaging values for each bin
%% value
homeData='/afs/cs/group/reconstruction3d/Data/Dataset_Oct21/scenario';
depthDirectory = strcat(homeData,num2str(dataNumber),'/');
imgDirectory = strcat(homeData,num2str(dataNumber),'/');
%Win05:
% system([ 'cat ' depthDirectory depthFilename ' | sed s/PTLASER/\/g | cat > ' depthDirectory depthFilename 'cleaned']);
% Above command is passed in Linux to clear the text, leaving only numbers, which can be directly
% imported into Matlab by using File > Import Data
    dirList = dir([depthDirectory '*.txt']);
    depthFileInitialname = dirList(1).name;
    
    system([ 'cat ' depthDirectory depthFileInitialname ' | sed /PTLASER/\d | cat > ' depthDirectory depthFileInitialname 'image1']);
    system([ 'cat ' depthDirectory depthFileInitialname 'image1' ' | sed s/IMAGE/\/g | cat > ' depthDirectory depthFileInitialname 'image2']);
    system([ 'cat ' depthDirectory depthFileInitialname 'image2' ' | sed s/data/\\n\data/g | cat > ' depthDirectory depthFileInitialname 'image3']);
    system([ 'cat ' depthDirectory depthFileInitialname 'image3' ' | sed /data/\d | cat > ' depthDirectory depthFileInitialname 'image']);
    system([ 'cat ' depthDirectory depthFileInitialname ' | sed s/PTLASER/\/g | cat > ' depthDirectory depthFileInitialname 'cleaned1']);
    system([ 'cat ' depthDirectory depthFileInitialname 'cleaned1' ' | sed /IMAGE/\d | cat > ' depthDirectory depthFileInitialname 'cleaned']);
    system([ 'rm ' depthDirectory depthFileInitialname 'cleaned1']);
    system([ 'rm ' depthDirectory depthFileInitialname 'image1']);
    system([ 'rm ' depthDirectory depthFileInitialname 'image2']);
    system([ 'rm ' depthDirectory depthFileInitialname 'image3']);
    
    depthFilename = strcat(depthFileInitialname,'cleaned');
    imageFileForAnglesName = strcat(depthFileInitialname,'image');
rawDepth = load([depthDirectory depthFilename]); %load([directory 'log' num2str(dataNumber) '.mat']);
imageAngList = load([depthDirectory imageFileForAnglesName]);
%% now we make the bins and store the rawDepth values in the appropriate
%% bins
panStartAngD  =   floor(rawDepth(1,2)-1);
panEndAngD    =   ceil(rawDepth(end,2)+1);
sweepStartAngD=   -90;
sweepEndAngD  =   90;
panRange = panEndAngD-panStartAngD;
sweepRange = sweepEndAngD-sweepStartAngD;
 
panBinResD  =   .125; %panRange/size(rawDepth,1);
sweepBinResD=   1; %sweepRange/rawDepth(1,4);
panBins     = ceil(panRange/panBinResD);  
sweepBins   = ceil(sweepRange/sweepBinResD);
depthBins = zeros(panBins,sweepBins);
countBuf = zeros(panBins,sweepBins);
[depthBins,countBuf,minDepth]=readLogFileForDepth(rawDepth,depthBins,countBuf, ...
    panEndAngD,panBinResD,sweepEndAngD,sweepBinResD);
%% trim depth buffer to remove zero columns
[depthBins,countBuf,panStartTrim,panEndTrim,sweepStartTrim,sweepEndTrim]=trimDepthBuf(depthBins,countBuf,panBins,sweepBins);
panBins=size(depthBins,1);
sweepBins=size(depthBins,2);
panStartAngD  =   floor(rawDepth(1,2))+(panStartTrim*panBinResD);%rawDepth(1,2);
panEndAngD    =   ceil(rawDepth(end,2))-(panEndTrim*panBinResD);%rawDepth(end,2);
sweepStartAngD=   -90+(sweepStartTrim*sweepBinResD);
sweepEndAngD  =   90-(sweepEndTrim*sweepBinResD);
panRange = panEndAngD-panStartAngD;
sweepRange = sweepEndAngD-sweepStartAngD;
%% we now check for validity of the buffer
[depthBins,countBuf]=checkIfValid(depthBins,countBuf,size(depthBins,1), ...
    size(depthBins,2),minDepth);
clear countBuf;
%% we now do the triangle making
%% for the time being we are skipping this part
%% now for every i,j in the depth map, we need to find the x_im and y_im,
%% based on the format that Min wants it in i.e. lets say that the image
%% lines between (0,0), (1,0), (1,1) and (0,1) going counter clockwise with
%% the bottom left point of the image being (0,0), and so the top left
%% point being (0,1)
% image2LaserPanAngleOffsetD = 0; %% check if this is 0 by trial and error
% panImageStartAngD = -180;
% panImageEndAngD = 180;
% sweepImageStartAngD = 67-90;
% sweepImageEndAngD = (180-68)-90;
% panBinIdxStart=floor((panEndAngD-panImageStartAngD)/panBinResD);
% panBinIdxEnd=floor((panEndAngD-panImageEndAngD)/panBinResD);
% sweepBinIdxStart=floor((sweepEndAngD-sweepImageStartAngD)/sweepBinResD);
% sweepBinIdxEnd=floor((sweepEndAngD-sweepImageEndAngD)/sweepBinResD);
        
laserPanAxisOffset = 0.2;
laserTiltAxisOffset = 0.155;
rigHeight = 1.5;
distCoeffs = [-0.068382333138247237 -0.54357925933026252 -0.0029516085790312376 ...
        -0.00024597914695977276];
D2R = (pi/180);
fx = 2400.2091651084;
fy = 2407.3312729885838;
Ox = 1110.7122391785729;%2272/2; %
Oy = 833.72104535435108;%1704/2; %
camIntMat = [ -fx, 0.0, Ox;
    0.0, -fy, Oy;
    0.0, 0.0, 1.0 ];
UNDISTORT = 1;
rollInit = pi;
pitchInit = -pi/2;
yawInit = 0.0;
camPanAxisOffset = 0.20;
camTiltAxisOffset = 0.155;
imagePanAngD = imageAngList(:,2);
imagePanAngD = (imagePanAngD<-180).*(imagePanAngD+360)+...
    (imagePanAngD>-180).*(imagePanAngD);
imagePanAngRndD = round(imageAngList(:,2));
imageRangeAngD = 21;
R = eulerAngles(rollInit,pitchInit,yawInit);
%% roll is along the x-axis, pitch is along the y-axis and yaw is along
%% the z axis
t = [camPanAxisOffset 0.0 camTiltAxisOffset];
camPos = [R t'; 0 0 0 1];
camExtMat1 = inv(camPos);
for imNum=1:length(imagePanAngD)
    image2LaserPanAngleOffsetD = 0; %% check if this is 0 by trial and error
    dirList = dir([imgDirectory strcat('*', ...
        num2str(abs(imagePanAngRndD(imNum))),'*.jpg')]);
    if (length(dirList)==0)
        problem_list=1
    end
	imageFilename = dirList(1).name;
    depthPixelFilename=regexprep(imageFilename,'img','depth_sph_corr');
    depthPixelFilename=regexprep(depthPixelFilename,'.jpg','.mat');
    depthPixelFilename=strcat('/afs/cs/group/reconstruction3d/scratch/LaserDepth/Dataset_Oct21/',depthPixelFilename);
    
    imagePanAngCurrentD=(round(imagePanAngD(imNum)/panBinResD))*panBinResD;
	panImageStartAngD = imagePanAngCurrentD-imageRangeAngD;
	panImageEndAngD = imagePanAngCurrentD+imageRangeAngD;
	sweepImageStartAngD = 63-90;
	sweepImageEndAngD = (180-63)-90;
    pixelCount = 1;
% 	panCount=1;
% 	imageFilenamePrev=' ';
	for panBinAngleLoopD=panImageEndAngD:-panBinResD:panImageStartAngD
        panBinAngleLoopDt=panBinAngleLoopD;
        if panBinAngleLoopD<-180 %% wrap around
            panBinAngleLoopDt=panBinAngleLoopD+360;
        end
        if panBinAngleLoopD>180 %% wrap around
            panBinAngleLoopDt=panBinAngleLoopD-360;
        end
        if (((panStartAngD<-180) && (panEndAngD>180))|| ...
                ((panBinAngleLoopD>=panStartAngD) && (panBinAngleLoopD<panEndAngD)))
            %% this condition is to handle log files which do not have the
            %% full pan angle from -180 to 180
        
            panBinIdx=floor((panEndAngD-panBinAngleLoopDt)/panBinResD);
    %         sweepCount=1;
            for sweepBinAngleLoopD=sweepImageStartAngD:sweepBinResD:sweepImageEndAngD
                sweepBinIdx=floor((sweepEndAngD-sweepBinAngleLoopD)/sweepBinResD);
                %% we now calculate the world coordinate point from depth buffer
                if ((panBinIdx>size(depthBins,1))||(panBinIdx<1)||...
                        (sweepBinIdx>size(depthBins,2))||(sweepBinIdx<1))
                    waithere = 1
                end           
                depthVal = depthBins(panBinIdx,sweepBinIdx);
                Wcoord=calcPointFromDepthBuf(depthVal,panBinIdx,sweepBinIdx, ...
                    panEndAngD,panBins,panRange,sweepEndAngD,sweepBins, ...
                    sweepRange,laserPanAxisOffset,laserTiltAxisOffset,rigHeight);
    % 			%% once we have the world coordinates, we get the color
    % 			%% corresponding to these world coordinates
    % 			%% which is the same as finding the x_im and y_im in the correct
    % 			%% image		
    % 			%% so the first step will be to find the right image
    %             currPanAngD = panEndAngD-(panBinIdx/panBins)*panRange;
    %             currPanAngD = (currPanAngD+image2LaserPanAngleOffsetD);
    %             if (currPanAngD>180)
    %                 currPanAngD = currPanAngD-360;
    %             elseif (currPanAngD<-180)
    %                 currPanAngD = currPanAngD+360;
    %             end
    %             [temp1,temp2]=min(abs(imagePanAngD-currPanAngD));
    %             rightImageNum = imagePanAngRndD(temp2);
        %         if (imageFilenamePrev~=imageFilename)
        %             A = imread([imgDirectory imageFilename]);
        %             imageFilenamePrev=imageFilename;            
        %         end
                %% the definition of roll, pitch and yaw needs to be checked with
                %% Kyle
                panAngR = imagePanAngD(imNum)*D2R;
                tiltAngR = imageAngList(imNum,3)*D2R;
                [pixel, camCoord, isCorrect] = getPixelCoords(Wcoord,panAngR,tiltAngR, ...
                    rigHeight,distCoeffs,camIntMat,UNDISTORT,camExtMat1);
                %% the storage format in the pixelreservoir is to store for each
                %% value of panBinIdx and sweepBinIdx, the imageNumber and the
                %% pixel values in that image
                if (isCorrect == true)                
                    pixelReservoir(pixelCount,:)=[camCoord' pixel depthVal];
                    pixelCount=pixelCount+1;                
                end
            end
        end
%             waithere=1;
%         panCount=panCount+1;
	end
    depthMapXRes=55;
	depthMapYRes=305;
%     [pseudoDepthMap] = sortPixel(pixelReservoir,depthMapXRes,depthMapYRes);
    save(depthPixelFilename, 'pixelReservoir');
%         save(depthPixelFilename, 'pseudoDepthMap', 'pixelReservoir');
%     pixelCount
% 	plotKyleData(pseudoDepthMap,depthMapXRes,depthMapYRes, ...
%         imgDirectory,imageFilename)
    
    clear pixelReservoir;
end
% displayScenario(depthDirectory,imgDirectory);
