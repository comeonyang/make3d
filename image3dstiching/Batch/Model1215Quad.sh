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
matlab -nojvm -nodisplay > Model1215Quad.out << EOF

% Run your MATLAB commands inlin
cd ..
Fdir ='/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/1215';
WlrName = '1215Quadtest_NewFarestDist'

Flag.NewInput = 0;
Flag.FlagDisp = 0;
Flag.PoseMatches = 0;
Flag.ReInference = 1;
Flag.PoseMatchStor = 0;
Flag.FindOcclu = 1;
Flag.NewModelList = 1; % if 0 used the prestoraged tri-angulation data
Flag.FlagFirstPair = Flag.NewModelList; 
Flag.FlagPreloadOccluDetect = 1;

PairList = [];
PairList = [ PairList; {'IMG_0128','IMG_0129'}];
%PairList = [ PairList; {'IMG_0128','IMG_0130'}];
%PairList = [ PairList; {'IMG_0128','IMG_0132'}];
PairList = [ PairList; {'IMG_0129','IMG_0130'}];
%PairList = [ PairList; {'IMG_0129','IMG_0132'}];
%PairList = [ PairList; {'IMG_0130','IMG_0132'}];
image3dstichingSetPath

tic;
main( Fdir, PairList, WlrName, 1, Flag);
toc;

cd Batch/;

% Exit MATLAB
exit
EOF

# Display the output
cat Model1215Quad.out
