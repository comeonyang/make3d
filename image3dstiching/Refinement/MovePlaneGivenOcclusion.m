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
function []=MovePlaneGivenOcclusion...
		(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, Img1, Img2, Wrlname, PostFixStrAfter,...
		 FaceSetPickedIND1, OccluedFaceSetIND1, FaceSetPickedIND2, OccluedFaceSetIND2, Mask1, Mask2)

% Pick the redundent plane and make the corresponding Sup to zero
	[model1 model2] = RedundentPlanePick(defaultPara, ImgInfo1, ImgInfo2, Pair, GlobalScale, ...
				FaceSetPickedIND1(Mask1), OccluedFaceSetIND1(Mask1), ...
				FaceSetPickedIND2(Mask2), OccluedFaceSetIND2(Mask2));

	% storage the new ImgInfo.Model.PlaneParaInfo.SupOri
	model = model1;
	save([defaultPara.Fdir '/data/' Img1 '/' Wrlname '_' Img1 '_' PostFixStrAfter '.mat'],'model');
	model = model2;
	save([defaultPara.Fdir '/data/' Img2 '/' Wrlname '_' Img2 '_' PostFixStrAfter '.mat'],'model');	
	
	return;
