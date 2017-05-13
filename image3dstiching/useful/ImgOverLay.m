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
function [Box]=ImgOverLay(Target, Field, FieldStitchPosition, TargetStitchPosition)

% []=ImgOverLay(Target, Field, FieldStitchPosition, TargetStitchPosition)
% This functino overlay the field and target image given at input
FieldStitchPosition = round( FieldStitchPosition);
TargetStitchPosition = round( TargetStitchPosition);

NumberOfChannel = size(Target,3);
TargetSize = size( Target);
FieldSize = size( Field);

BoxBottomRight = max( [ size(Field); FieldStitchPosition + ( TargetSize - TargetStitchPosition)], [], 1);
BoxTopLeft = min( [ ones(1,2); FieldStitchPosition - TargetStitchPosition + 1], [], 1);
BoxSize = BoxBottomRight - BoxTopLeft + 1;

FieldTopLeftPosition = ones(1,2) - BoxTopLeft + 1;
TargetTopLeftPosition = ( FieldStitchPosition - TargetStitchPosition + 1) -BoxTopLeft + 1;

% initialize the Box holding Target and Field
BoxTarget = repmat( zeros( BoxSize), [1 1 NumberOfChannel]);
BoxField = repmat( zeros( BoxSize), [1 1 NumberOfChannel]);
Box = zeros(size(BoxField));

% change the order of dimension
BoxTarget = permute(BoxTarget, [3 1 2]);
BoxField = permute(BoxField, [3 1 2]);
Box = permute(Box, [3 1 2]);
Target = permute(Target, [3 1 2]);
Field = permute(Field, [3 1 2]);

% enter info of Taget and Field in specific position in BoxTarget and BoxField
BoxField(:,FieldTopLeftPosition(1):( FieldTopLeftPosition(1)+FieldSize(1)-1), FieldTopLeftPosition(2):(FieldTopLeftPosition(2)+FieldSize(2)-1)) = Field;
BoxTarget(:,TargetTopLeftPosition(1):( TargetTopLeftPosition(1)+TargetSize(1)-1), TargetTopLeftPosition(2):(TargetTopLeftPosition(2)+TargetSize(2)-1)) = Target;

% Merging this two images
MaskTarget = BoxTarget(1,:,:) ~= 0 & BoxField(1,:,:) ==0;
MaskField = BoxTarget(1,:,:) == 0 & BoxField(1,:,:) ~=0;
for i = 1:NumberOfChannel
    Box(i,MaskTarget) = BoxTarget(i,MaskTarget);
    Box(i,MaskField) = BoxField(i,MaskField);
    Box(i,~MaskTarget & ~MaskField) = (BoxField(i,~MaskTarget & ~MaskField) + BoxTarget(i,~MaskTarget & ~MaskField))/2;
end    

Box = permute(Box, [2 3 1]);
return;
