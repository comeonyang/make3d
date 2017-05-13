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


% This is a script to use to build tedous data building


I1=imreadbw([Fdir '/pgm/' PairImg '.pgm']); % function from sift
I2=imreadbw([Fdir '/pgm/' NewImg '.pgm']); % function from sift

[f1, f2] = readSurf(PairImg, NewImg, Fdir, 1);

hit = 0;
for i = 1:length(ImgInfo)
	if strcmp(ImgInfo(i).ExifInfo.name, [PairImg '.jpg'] )
		ImgInfo1 = ImgInfo(i);
		PairImgIndex = i;
		hit = hit + 1;
	elseif strcmp(ImgInfo(i).ExifInfo.name, [NewImg '.jpg'] )
		ImgInfo2 = ImgInfo(i);
		NewImgIndex = i;
		hit = hit + 1;
	elseif hit >=2
		break;
	end
end

%depth
[DV DH] = size(Depth2);
[IV IH] = size(I2);
Scale = [[DV DH];[IV IH]];
[IND2]=ReScaleInd2Sub(Scale,f2); 
D2 = Depth2(IND2);

% load pre_data in ImgAddedList
load([Fdir '/data/' PairImg '_Data.mat']);
Depth1 = Depth;
Sup1 = Sup;
[DV DH] = size(Depth1);
[IV IH] = size(I1);
Scale = [[DV DH];[IV IH]];
[IND1]=ReScaleInd2Sub(Scale,f1);
D1 = Depth1(IND1);


