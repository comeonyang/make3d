# *  This code was used in the following articles:
# *  [1] Learning 3-D Scene Structure from a Single Still Image, 
# *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
# *      In ICCV workshop on 3D Representation for Recognition (3dRR-07), 2007.
# *      (best paper)
# *  [2] 3-D Reconstruction from Sparse Views using Monocular Vision, 
# *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
# *      In ICCV workshop on Virtual Representations and Modeling 
# *      of Large-scale environments (VRML), 2007. 
# *  [3] 3-D Depth Reconstruction from a Single Still Image, 
# *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
# *      International Journal of Computer Vision (IJCV), Aug 2007. 
# *  [6] Learning Depth from Single Monocular Images, 
# *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
# *      In Neural Information Processing Systems (NIPS) 18, 2005.
# *
# *  These articles are available at:
# *  http://make3d.stanford.edu/publications
# * 
# *  We request that you cite the papers [1], [3] and [6] in any of
# *  your reports that uses this code. 
# *  Further, if you use the code in image3dstiching/ (multiple image version),
# *  then please cite [2].
# *  
# *  If you use the code in third_party/, then PLEASE CITE and follow the
# *  LICENSE OF THE CORRESPONDING THIRD PARTY CODE.
# *
# *  Finally, this code is for non-commercial use only.  For further 
# *  information and to obtain a copy of the license, see 
# *
# *  http://make3d.stanford.edu/publications/code
# *
# *  Also, the software distributed under the License is distributed on an 
# * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# *  express or implied.   See the License for the specific language governing 
# *  permissions and limitations under the License.
# *
# */

#!/bin/bash

# Change to the submission directory
cd $PBS_O_WORKDIR

# Run the m-file
matlab -nojvm -nodisplay > CleanMatches.out << EOF

% Run your MATLAB commands inlin
cd ..
if false
k = 1;
PairList = [];
for i = 1:5%length(ImgInfo)
	for j = (i+1):length(ImgInfo)
		PairList{k,1} = strrep(ImgInfo(i).ExifInfo.name, '.jpg', '');
		PairList{k,2} = strrep(ImgInfo(j).ExifInfo.name, '.jpg', '');
		k = k + 1;
	end
end
else
	if false
	PairList = [ PairList; {'IMG_0664','IMG_0668'}];
	PairList = [ PairList; {'IMG_0668','IMG_0674'}];
	PairList = [ PairList; {'IMG_0674','IMG_0676'}];
	PairList = [ PairList; {'IMG_0676','IMG_0677'}];
	PairList = [ PairList; {'IMG_0677','IMG_0678'}];

	% inlcude all posible pairs
	PairList = [ PairList; {'IMG_0663','IMG_0665'}];
	PairList = [ PairList; {'IMG_0663','IMG_0667'}];
	PairList = [ PairList; {'IMG_0663','IMG_0669'}];
	PairList = [ PairList; {'IMG_0663','IMG_0670'}];
	PairList = [ PairList; {'IMG_0663','IMG_0671'}];

	PairList = [ PairList; {'IMG_0664','IMG_0665'}];
	PairList = [ PairList; {'IMG_0664','IMG_0666'}];
	PairList = [ PairList; {'IMG_0664','IMG_0674'}];
	PairList = [ PairList; {'IMG_0664','IMG_0675'}];
	PairList = [ PairList; {'IMG_0664','IMG_0678'}];

	PairList = [ PairList; {'IMG_0665','IMG_0667'}];
	PairList = [ PairList; {'IMG_0665','IMG_0667'}];
	PairList = [ PairList; {'IMG_0665','IMG_0671'}];
	PairList = [ PairList; {'IMG_0665','IMG_0677'}];
	PairList = [ PairList; {'IMG_0665','IMG_0679'}];

	PairList = [ PairList; {'IMG_0666','IMG_0667'}];
	PairList = [ PairList; {'IMG_0666','IMG_0668'}];
	PairList = [ PairList; {'IMG_0666','IMG_0671'}];
	PairList = [ PairList; {'IMG_0666','IMG_0677'}];
	PairList = [ PairList; {'IMG_0666','IMG_0674'}];
	end
%-----------------
	if false
	PairList = [ PairList; {'IMG_0667','IMG_0670'}];
	PairList = [ PairList; {'IMG_0667','IMG_0673'}];
	PairList = [ PairList; {'IMG_0667','IMG_0674'}];
	PairList = [ PairList; {'IMG_0667','IMG_0677'}];
	PairList = [ PairList; {'IMG_0667','IMG_0678'}];

	PairList = [ PairList; {'IMG_0668','IMG_0670'}];
	PairList = [ PairList; {'IMG_0668','IMG_0673'}];
	PairList = [ PairList; {'IMG_0668','IMG_0674'}];
	PairList = [ PairList; {'IMG_0668','IMG_0677'}];
	PairList = [ PairList; {'IMG_0668','IMG_0678'}];

	PairList = [ PairList; {'IMG_0669','IMG_0670'}];
	PairList = [ PairList; {'IMG_0669','IMG_0671'}];
	PairList = [ PairList; {'IMG_0669','IMG_0673'}];

	PairList = [ PairList; {'IMG_0670','IMG_0673'}];
	PairList = [ PairList; {'IMG_0670','IMG_0674'}];
	PairList = [ PairList; {'IMG_0670','IMG_0683'}];
	PairList = [ PairList; {'IMG_0670','IMG_0684'}];

	PairList = [ PairList; {'IMG_0671','IMG_0673'}];
	PairList = [ PairList; {'IMG_0671','IMG_0674'}];
	PairList = [ PairList; {'IMG_0671','IMG_0683'}];
	PairList = [ PairList; {'IMG_0671','IMG_0684'}];

	PairList = [ PairList; {'IMG_0673','IMG_0674'}];
	PairList = [ PairList; {'IMG_0673','IMG_0683'}];
	PairList = [ PairList; {'IMG_0673','IMG_0684'}];
	end
%----------------------
	if true
	PairList = [ PairList; {'IMG_0674','IMG_0676'}];
	PairList = [ PairList; {'IMG_0674','IMG_0677'}];
	PairList = [ PairList; {'IMG_0674','IMG_0678'}];
	PairList = [ PairList; {'IMG_0674','IMG_0679'}];

	PairList = [ PairList; {'IMG_0675','IMG_0678'}];
	PairList = [ PairList; {'IMG_0675','IMG_0677'}];
	PairList = [ PairList; {'IMG_0675','IMG_0678'}];
	PairList = [ PairList; {'IMG_0675','IMG_0679'}];

	PairList = [ PairList; {'IMG_0676','IMG_0677'}];

	PairList = [ PairList; {'IMG_0677','IMG_0678'}];
	PairList = [ PairList; {'IMG_0677','IMG_0679'}];

	PairList = [ PairList; {'IMG_0678','IMG_0679'}];
	PairList = [ PairList; {'IMG_0678','IMG_0680'}];
	PairList = [ PairList; {'IMG_0678','IMG_0681'}];

	PairList = [ PairList; {'IMG_0679','IMG_0680'}];
	PairList = [ PairList; {'IMG_0679','IMG_0681'}];
	PairList = [ PairList; {'IMG_0679','IMG_0682'}];

	PairList = [ PairList; {'IMG_0683','IMG_0684'}];
	end
end
EvaluateMatchesTracks
cd Batch/;

% Exit MATLAB
exit
EOF

# Display the output
cat CleanMatches.out
