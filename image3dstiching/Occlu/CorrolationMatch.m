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
function [Matches CoeffM Inliers]=CorrolationMatch( defaultPara, Pair, ITarget, IField, PointPix, POriReprojM, FieldOccluPix, MinMaxScaleFactor, Target2NumberOfPatches)

% Find corrolation match along the epipolar line
% by calling FindTarget.m for each pair of search 

% initializing Parameters
% MatchTHre = 0.3; % Magic number change later
if nargin <9
	%Target2NumberOfPatches = [ 2.5 5];
	Target2NumberOfPatches = [ 1 1]; % used smaller target size to speed up
end
ScaleFactorStep = 0.1;
TargetReduceRatio = 0.95;
[Iy Ix] = size(ITarget);
[IyField IxField] = size(IField);
TargetSize = round([Iy/55*Target2NumberOfPatches(1) Ix/305*Target2NumberOfPatches(2)]);% prevous used 5 and 10
NumTarget = size(PointPix,2);

% Process the Target_tmp region
Target_tmp_top = round( PointPix(2,:) - TargetSize(1) );
Target_tmp_bottom = round( PointPix(2,:) + TargetSize(1) );
Target_tmp_left = round( PointPix(1,:) - TargetSize(2) );
Target_tmp_right = round( PointPix(1,:) + TargetSize(2) );

% detect region exceed the boundary of images and shift the region accordingly
%	Yet make sure the PointPix is as close to the center of the Target_tmp as possible
% (Shift the Target_tmp vertically)
Mask =  Target_tmp_top < 1;
shift_vertical(Mask) = (1 - Target_tmp_top(Mask));
Target_tmp_top(Mask) = Target_tmp_top(Mask) + shift_vertical(Mask);
Target_tmp_bottom(Mask) = Target_tmp_bottom(Mask) + shift_vertical(Mask);
Mask = Target_tmp_bottom > Iy;
shift_vertical(Mask) = ( Iy - Target_tmp_bottom(Mask));
Target_tmp_top(Mask) = Target_tmp_top(Mask) + shift_vertical(Mask);
Target_tmp_bottom(Mask) = Target_tmp_bottom(Mask) + shift_vertical(Mask);
% (Shift the Target_tmp horizontal)
Mask = Target_tmp_left < 1;
shift_horizontal(Mask) = (1 - Target_tmp_left(Mask));
Target_tmp_left(Mask) = Target_tmp_left(Mask) + shift_horizontal(Mask);
Target_tmp_right(Mask) = Target_tmp_right(Mask) + shift_horizontal(Mask);
Mask = Target_tmp_right > Ix;
shift_horizontal(Mask) = (Ix - Target_tmp_right(Mask));
Target_tmp_left(Mask) = Target_tmp_left(Mask) + shift_horizontal(Mask);
Target_tmp_right(Mask) = Target_tmp_right(Mask) + shift_horizontal(Mask);

% Define the TopMargin, BottomMargin, LeftMargin, and RightMargin of Target_tmp
PointPixTopMargin = PointPix(2,:) - Target_tmp_top; % All Margin are positive
PointPixBottomMargin =  Target_tmp_bottom - PointPix(2,:);
PointPixLeftMargin =  PointPix(1,:) - Target_tmp_left; 
PointPixRightMargin =  Target_tmp_right - PointPix(1,:); 

% initializing variables
Matches = zeros(4,NumTarget);
%CoeffM = zeros(1,NumTarget);
CoeffM = zeros(2,NumTarget);

% 0) determine all possible Scale and region
ScaleFactors = MinMaxScaleFactor(1):ScaleFactorStep:MinMaxScaleFactor(2);
	
% 1) Field that fully contain target_tmp	
[ Rc ConS RoughConS tempConS] = EndPoint2BoxConS(defaultPara, Ix, Iy, POriReprojM, FieldOccluPix, 1);
Field_top = round( RoughConS(3,:) - PointPixTopMargin*MinMaxScaleFactor(2));
Field_left = round( RoughConS(1,:) - PointPixLeftMargin*MinMaxScaleFactor(2));
Field_bottom =  round( RoughConS(4,:)  + PointPixBottomMargin*MinMaxScaleFactor(2));
Field_right = round( RoughConS(2,:) + PointPixRightMargin*MinMaxScaleFactor(2));

tic
for i= 1:NumTarget
    if rem(i,100) == 0
        disp([ num2str( floor(i/100)) 'Target Tested in ' num2str(NumTarget)])
    end    
	% Pick the Proper target_tmp
	ITarget_tmp = ITarget( Target_tmp_top(i):Target_tmp_bottom(i), Target_tmp_left(i):Target_tmp_right(i) );

	% Define the Point to pick inside the Target Patch
	[ PickRowInTarget_tmp PickColInTarget_tmp] = PickAnchorPoint(ITarget_tmp);
	PickPointTarget(:,i) = [ PickRowInTarget_tmp + Target_tmp_top(i) - 1; PickColInTarget_tmp+ Target_tmp_left(i) - 1];

	% prepare to run FindTarget

	% Stop finfing Matches if the IField is crop by the Field image region, since it might cause Target not totally inside the Field
	if Field_top(i) < 1 || Field_left(i) < 1 || Field_bottom(i) > IyField || Field_right(i) > IxField
        if defaultPara.Flag.FlagCorrRefinement
            Field_top(i)  = min( max( Field_top(i), 1), IyField);
            Field_bottom(i)  = min( max( Field_bottom(i), 1), IyField);
            Field_left(i)  = min( max( Field_left(i), 1), IxField);
            Field_right(i)  = min( max( Field_right(i), 1), IxField);
        else    
       		continue;
        end    
	end

	% Pick the Field_tmp
	IField_tmp = IField(Field_top(i):Field_bottom(i), Field_left(i):Field_right(i));

	% 2) Epipolar geometry
        T1_hat = [[0 -Pair.T(3) Pair.T(2)];...
            [Pair.T(3) 0 -Pair.T(1)];...
            [-Pair.T(2) Pair.T(1) 0]];
        F = inv(defaultPara.InrinsicK2)'*T1_hat*Pair.R*inv(defaultPara.InrinsicK1);
	EpipoalLine = F*[ PointPix; ones(1,size(PointPix,2))];
	NormalizedEpipoalLine = EpipoalLine./repmat( sqrt(sum( EpipoalLine(1:2,:).^2,1)),3,1);
% if i == 10 || i ==34 || i ==45
%    pause 
% end    
	% Run Find Target
	[ X, Y, Scale, Coeff, Coeffs ] = FindTarget( defaultPara, IField_tmp, ITarget_tmp, Field_top(i), Field_left(i), PickRowInTarget_tmp, PickColInTarget_tmp, NormalizedEpipoalLine(:,i), ...
                [Iy Ix], ScaleFactors, 1, 1);% Min disable EpipolarLineFlag
	X = X+Field_left(i);
	Y = Y+Field_top(i);
	CoeffM(:,i) = Coeff;
	Matches(:,i) = [ PickPointTarget([2 1],i); X; Y];
%      if i == 10 || i ==34 || i ==45
% 	target_tmp = imresize(ITarget_tmp,Scale,'bicubic');
% 	TargetPickPoint = [PickRowInTarget_tmp PickColInTarget_tmp]*Scale;
% 	Box= ImgOverLay(target_tmp, IField, [Y X], TargetPickPoint);
% 	figure(1); imshow(Box);
% 	hold on; scatter(X,Y,'r');hold off;
% 	figure(2); imshow(ITarget);
% 	hold on; scatter( PickPointTarget(2,i), PickPointTarget(1,i),'r');hold off;
%     pause;
%      end
end

% Clean the redundent matches so that each PointPix match a unique points in the Field image //Min add July 10th================
[Inliers] = CleanMatch(Matches, CoeffM(1,:));
[InliersReverse] = CleanMatch(Matches(:,Inliers), CoeffM(1,Inliers));
Inliers = Inliers(InliersReverse);
% ===============================================================================================================================

TimeForCorrolationMatch = toc;
if defaultPara.Flag.FlagRefinementDisp
	disp([ 'End Of CorrolationMatch, it takes ' num2str(TimeForCorrolationMatch) ' seconds']);
end

return;
