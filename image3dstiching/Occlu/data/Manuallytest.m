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
AnchorXTarget = 876;
AnchorYTarget = 591;
AnchorXField = 1189;
AnchorYField = 780;
TargetSize = [100 150];
FieldRatio = 3;
TargetScale = 1.5;
GuessTargetPosition = [97 78];

ITarget_tmp = GrayI606(AnchorXTarget:(AnchorXTarget+TargetSize(1)),AnchorYTarget:(AnchorYTarget+TargetSize(2)));
figure; imshow(ITarget_tmp)

IField_tmp = GrayI610(AnchorYField:(AnchorYField+TargetSize(1)*FieldRatio),AnchorXField:(AnchorXField+TargetSize(2)*FieldRatio));
figure; imshow(IField_tmp)

% example overlay to decide the scale
Target2Field = zeros(size(IField_tmp));
Target2Field(GuessTargetPosition(1):(GuessTargetPosition(1)+size(ITarget_tmp,1)-1),GuessTargetPosition(2):(GuessTargetPosition(2)+size(ITarget_tmp,2)-1)) = ITarget_tmp;
Target2Field = uint8(Target2Field);
figure; imshow(Target2Field)

% Recale and do Cross_Corrolation match
ITarget_tmp_resized = imresize(ITarget_tmp, TargetScale, 'bicubic');
IField_tmp = GrayI610(AnchorYField:(AnchorYField+size(ITarget_tmp_resized,1)*FieldRatio),AnchorXField:(AnchorXField+size(ITarget_tmp_resized,2)*FieldRatio));
figure; imshow(IField_tmp)

res             = normxcorr2( ITarget_tmp_resized, IField_tmp );
res_top = size(ITarget_tmp_resized,1);
res_left = size(ITarget_tmp_resized,2);
res_modified = res(res_top:(end-res_top-1),res_left:(end-res_left-1));
figure; imagesc(abs(res_modified))

% find max position
[VCol ICol] = max(res_modified);
[V I] = max(VCol)
ICol(I)

Target2Field = zeros(size(IField_tmp));
Target2Field(ICol(I):(ICol(I)+size(ITarget_tmp_resized,1)-1),I:(I+size(ITarget_tmp_resized,2)-1)) = ITarget_tmp_resized;
Target2Field = uint8(Target2Field);
Ioverlay = cat(3, IField_tmp,Target2Field,zeros(size(IField_tmp)));
figure; imagesc(Ioverlay)
