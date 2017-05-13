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
ImgTarget = ImgInfo(1).ExifInfo.IDName;
ITarget = imreadbw([defaultPara.Fdir '/pgm/' ImgTarget '.pgm']);
ImgField = ImgInfo(2).ExifInfo.IDName;
IField = imreadbw([defaultPara.Fdir '/pgm/' ImgField '.pgm']);
TTarget_Field_hat = [[0 -Pair.T(3) Pair.T(2)];...
            [Pair.T(3) 0 -Pair.T(1)];...
            [-Pair.T(2) Pair.T(1) 0]];
FTarget_Field = inv(defaultPara.InrinsicK2)'*TTarget_Field_hat*Pair.R*inv(defaultPara.InrinsicK1);
MaxD1 = ones(1,length(TargetPointPix3_9));
MinD1 = ones(1,length(TargetPointPix3_9));
MaxD2 = ones(1,length(TargetPointPix9_3));
MinD2 = ones(1,length(TargetPointPix9_3));

dispMatchSearchRegin(ITarget, IField, [TargetPointPix3_9; ones(1,length(TargetPointPix3_9))], [TargetPointPix9_3; ones(1,length(TargetPointPix9_3))], Region9_3, Region3_9,...
            FTarget_Field, ...
            POriReprojM3_9, MaxD1, PoccluM3_9, MinD1, ...
            POriReprojM9_3, MaxD2, PoccluM9_3, MinD2, ...
            0, 'Stacking', 'v', 'Interactive', 0);
