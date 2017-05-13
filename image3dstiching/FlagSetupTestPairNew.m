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
Flag.NewInput = 0;
Flag.FlagDisp = 0;
Flag.PoseMatches = 1;
Flag.ReInference = 1;
Flag.PoseMatchStor = 1;
Flag.FindOcclu = 1;
Flag.NewModelList = 1; % if 0 used the prestoraged tri-angulation data
Flag.FlagFirstPair = Flag.NewModelList; % should equal Flag.NewModelList, since when newModelList is 1 then the first image is the FirstPair
Flag.FlagStorageBeforePairNew = 1;
Flag.FlagStorageAfterPairNew = 1;
% refinement in Metricon.m
Flag.Refinement = 1;
Flag.RefinementInference = 1;
Flag.FlagImgInfoLoadPreStorage = 4;
Flag.LoadStoragedSurfOccluMatches = 0;
Flag.LoadCorrMatches = 0;
Flag.StorageDataBeforeRansac = 1;

% Flag for Refinement step
Flag.FlagPreloadOccluDetectMatches = 0;
Flag.FlagPreloadOccluDetectRay = 0;
Flag.loadStoragedPairNew = 0;
Flag.FlagOccluDetect = 0;
Flag.FlagRefinementReInference = 1;
Flag.FlagRefinementDisp = 1;
Flag.FlagResetModeltoNomoModel = 1;
Flag.FlagRecipicalOccluDetection = 1;
Flag.FlagCorrRefinement = 0; % Not been used for a while July25th
Flag.UseNormXCorr2 = 2;
