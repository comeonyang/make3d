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

image3dstichingSetPath
Fdir = '/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/COM4_070525_001900_gates/data'
defaultPara.Fdir = '/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/COM4_070525_001900_gates'
load([ Fdir '/ImgInfo.mat']);
[ MatchesUnified]= MergeMatchesData(Fdir,{'MatchesType', 'MatchesType2','MatchesType3'});
% ----------- modified the match -
BadPairList = [];
BadPairList = [BadPairList; {'IMG_0666','IMG_0674'}];
BadPairList = [BadPairList; {'IMG_0670','IMG_0683'}];
BadPairList = [BadPairList; {'IMG_0670','IMG_0684'}];
BadPairList = [BadPairList; {'IMG_0671','IMG_0674'}];
BadPairList = [BadPairList; {'IMG_0671','IMG_0683'}];
BadPairList = [BadPairList; {'IMG_0671','IMG_0684'}];
BadPairList = [BadPairList; {'IMG_0673','IMG_0674'}];
BadPairList = [BadPairList; {'IMG_0673','IMG_0683'}];
BadPairList = [BadPairList; {'IMG_0674','IMG_0679'}];

GoodPairList = [];
GoodPairList = [ GoodPairList; {'IMG_0664','IMG_0668'}];
GoodPairList = [ GoodPairList; {'IMG_0668','IMG_0674'}];
GoodPairList = [ GoodPairList; {'IMG_0674','IMG_0676'}];
GoodPairList = [ GoodPairList; {'IMG_0676','IMG_0677'}];
GoodPairList = [ GoodPairList; {'IMG_0677','IMG_0678'}];
[MatchesUnified] = ModifyMatch(ImgInfo, MatchesUnified, GoodPairList, 1);
% --------------------------------
[ track]=TrackBuilding(MatchesUnified);
[ CleanedMatches]=Track2Match(track)

NuPair = size(GoodPairList,1);
for     i = 1:NuPair
	[ImgAIndex] = ImgInfoIndexFromName(ImgInfo, GoodPairList{i,1})
        [ImgBIndex] = ImgInfoIndexFromName(ImgInfo, GoodPairList{i,2})
	Img1 = strrep(ImgInfo(ImgAIndex).ExifInfo.name,'.jpg','');
	Img2 = strrep(ImgInfo(ImgBIndex).ExifInfo.name,'.jpg','');
	matches = [CleanedMatches(ImgAIndex, ImgBIndex).Index; CleanedMatches(ImgBIndex, ImgAIndex).Index ];
	matches_ori = [MatchesUnified(ImgAIndex, ImgBIndex).Index; MatchesUnified(ImgBIndex, ImgAIndex).Index ];

	I1=imreadbw([defaultPara.Fdir '/pgm/' Img1 '.pgm']); % function from sift
	I2=imreadbw([defaultPara.Fdir '/pgm/' Img2 '.pgm']); % function from sift
	[f1] = readSurf(Img1, defaultPara.Fdir, 'Dense'); % original features
	[f2] = readSurf(Img2, defaultPara.Fdir, 'Dense'); % original features

	figure(11);  plotmatches(I1,I2,f1, f2, matches, 'Stacking','v','Interactive', 2);
	saveas(11,[ Fdir '/cleanMatch/' Img1 '_' Img2 'CleanedMatch'],'jpg');
	figure(12);  plotmatches(I1,I2,f1, f2, matches_ori([2 1],:), 'Stacking','v','Interactive', 2);
	saveas(12,[ Fdir '/cleanMatch/' Img1 '_' Img2 'OriMatch'],'jpg');
end

