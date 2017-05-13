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
function makeDepthImages(dataNumber, displayFlag)

if nargin < 2
    displayFlag = false;
end
if nargin < 1
    dataNumber = 1;
end
centerAngle = -171.7 + dataNumber*31;

%Win05: Recalculate these numbers by looking at images and depthmaps 
laserAngleScan = 19;
topAngle = 63; %% earlier 63
botAngle = 63; %% if we make it 64, the left side of the image is more aligned, however if we make it 62, 
%% the right siede of the image is more aligned, so we are choosing it to be 63
leftAngle = 10;
rightAngle = 7;
horThetaStep = 0.125; %.5

% close all;
% depthDirectory = './dataset1/';
% imgDirectory = './dataset1/';
% 
% if uniqueCode > 0
%     dirList = dir([depthDirectory '*.txt']);
%     depthFilename = dirList(uniqueCode).name;

homeData='/afs/cs/group/reconstruction3d/Data/Dataset_July7/scenario';
depthDirectory = strcat(homeData,num2str(dataNumber),'/');
imgDirectory = strcat(homeData,num2str(dataNumber),'/');
%Win05:
% system([ 'cat ' depthDirectory depthFilename ' | sed s/PTLASER/\/g | cat > ' depthDirectory depthFilename 'cleaned']);
% Above command is passed in Linux to clear the text, leaving only numbers, which can be directly
% imported into Matlab by using File > Import Data

    dirList = dir([depthDirectory '*.txt']);
    depthFilename = dirList(1).name;
%     use to remove PTLASER text from files
%     system([ 'cat ' depthDirectory depthFilename ' | sed s/PTLASER/\/g | cat > ' depthDirectory depthFilename 'cleaned']);
    depthFilename = strcat(depthFilename,'cleaned');
%     dirList = dir([imgDirectory '*right*.jpg']);
    dirList = dir([imgDirectory '*.jpg']);
    %dirList = dir([imgDirectory 'packard1_10-28-05right*.jpg']);
    imNumber = -dataNumber+6;
    if imNumber <=0
        imNumber = imNumber+13;
        if imNumber == 7
            imNumber = 6;
        end
    end
    imageFilename = dirList(imNumber).name;
    %% 5 is 1 going down and coming back in a circle from 13, 7 is not
    %% there, but 6 is. So 1-5, 2-4, 3-3, 4-2, 5-1, 6-13, 7-12, 8-11, 9-10,
    %% , 10-9, 11-8, 12-6 skipping 7

logDepth = load([depthDirectory depthFilename]); %load([directory 'log' num2str(dataNumber) '.mat']);
A = imread([imgDirectory imageFilename]); %A = imread([directory 'img-' num2str(dataNumber) '.jpg']);
x = (centerAngle-laserAngleScan):horThetaStep:(centerAngle+laserAngleScan);

sizeHor = round( (-min(logDepth(:,2)) + max(logDepth(:,2) ) )/ horThetaStep );
quantizedLogDepth = zeros( length(x), size(logDepth,2)-4 );
angleIndex = 0;

for xc = (centerAngle-laserAngleScan):horThetaStep:(centerAngle+laserAngleScan)
    if (xc>max(logDepth(:,2)))
        x=(xc-180)-180;%212;
%         if (x<min(logDepth(:,2)))
%             x=max(logDepth(:,2))+x-min(logDepth(:,2));
%         end
    else
        x=xc;
    end
    [tmp ind] = min( abs( logDepth(:,2) - x) );
	%Win05: I remember 2nd column contained the time stamp, therefore it basically matches 
	% the timestamp
    %IMPROVEMENT: instead of nearest one, do AVERAGING
    angleIndex = angleIndex + 1;
    quantizedLogDepth(angleIndex,:) = logDepth(ind, 5:size(logDepth,2) );
	%Win05: 5 because first 4 numbers are not useful
end

clippedDepth = quantizedLogDepth( size(quantizedLogDepth,1):-1:1, ...
        (size(quantizedLogDepth,2)-topAngle):-1:botAngle )';

min(clippedDepth(:));
depthMap = (clippedDepth);    
%depthImg = clippedDepth .^ (1/3);    
% displayDepthMaps(depthMap);
% figure;
% imshow(A);

% if displayFlag
%     %subplot(1,2,1);
%     image(A);
%     axis square;
%     %subplot(1,2,2); 
%     figure,
%     imagesc( depthMap( :, leftAngle:(size(depthMap,2)-rightAngle) ) .^ (1/3) );
%     axis square;
% end
% 
% % image A = size(500,500,3 ), dephmap D = size(25,25), D_high = imresize(D, 500, 500, 'nearest')
% % B = rgb2ycbcr(A); B(:,:,2) = scaling factor * D_high, % set color channel tp depthmap.
% % imagesc(A), axis equal; imagesc(B), axis equal; displaydepthMap(D);
% % subplot(2,2,..)
% 
% 
% % size(A)
% % size(depthMap)
% % D_high = imresize(depthMap, [size(A,1) size(A,2)], 'nearest');
% % minDistance = min(min(D_high));
% % maxDistance = max(max(D_high));
% % B = rgb2hsv(A); 
% % scaling_factor = 1/maxDistance;
% % B(:,:,1) = scaling_factor*D_high;
% % % B(:,:,2) = scaling_factor*D_high; % The scaling factor should be chosen such that values are in between
% % %% 0 and 255? Here we set color channel to depthmap.
% % % 
% % figure;
% % subplot(1,2,1);
% % imagesc(A), axis equal; 
% % 
% % %% copying from displayDepthMap
% % 
% % jetMap = jet;
% % % close;
% % jetMap = 1-jetMap(end:-1:1,:);
% % jetMap = jetMap(end:-1:1,:);
% % jetMap = imresize(jetMap, [256 3], 'bilinear');
% % displayNorm = size(jetMap,1);
% % warning off;
% % D_high = uint8( displayNorm* (D_high - minDistance) / (maxDistance - minDistance) );
% % warning on;
% % subplot(1,2,2);
% % imagesc( D_high ); axis equal; colormap( jetMap );    axis off;   axis tight;
% % 
% % figure; %subplot(2,2,2);
% % imagesc(B), axis equal; 
% 
% size(A)
% size(depthMap)
% D_high = imresize(depthMap, [size(A,1) size(A,2)], 'nearest');
% % B = rgb2hsv(A); 
% % B = rgb2ycbcr(A);
% % scaling_factor = 3;
% % B(:,:,2) = scaling_factor*D_high; % The scaling factor should be chosen such that values are in between
% %% 0 and 255? Here we set color channel to depthmap.
% % 
% % figure;
% % subplot(1,2,1);
% % imagesc(A), axis equal; 
% 
% %% copying from displayDepthMap
% 
% % minDistance = min(min(D_high));
% % maxDistance = max(max(D_high));
% % jetMap = jet;
% % close;
% % jetMap = 1-jetMap(end:-1:1,:);
% % jetMap = jetMap(end:-1:1,:);
% % jetMap = imresize(jetMap, [256 3], 'bilinear');
% % displayNorm = size(jetMap,1);
% % warning off;
% % D_high = uint8( displayNorm* (D_high - minDistance) / (maxDistance - minDistance) );
% % warning on;
% % subplot(1,2,2);
% % imagesc( D_high ); axis equal; colormap( jetMap );    axis off;   axis tight;
% 
% % figure; %subplot(2,2,2);
% % imagesc(B), axis equal; 
% 
% % if dataNumber == 5
% %     waithere = 1;
% % end
% % displaydepthMap(D);
% % 
% % subplot(2,2,..);

directory = '/afs/cs/group/reconstruction3d/scratch/Min/rawlaserdata/';
save([directory strrep(strrep(imageFilename,'img','depth'),'.jpg','.mat')], 'depthMap' );
% save([directory strrep(strrep(imageFilename,'img','depth_high_res'),'.jpg','.mat')], 'D_high' );
% save([directory 'calculatedDepthData-' num2str(dataNumber) '.mat'], 'depthMap' );
% imwrite( depthMap/ max(depthMap(:) ), [directory 'calculatedDepthImgSet2-' num2str(dataNumber) '.jpg']);
