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
matlab -nojvm -nodisplay > Model.out << EOF

% Run your MATLAB commands inlin
cd ..
image3dstichingSetPath
Fdir ='/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/ImgWithLaserData';
WlrName = 'TestLaserImg'
PairList = [];

%        PairList = [ PairList; {'img-manmade2-p-169t0','img-manmade6-p-107t0'}];
       % PairList = [ PairList; {'img-manmade2-p-169t0','img-manmade7-p-343t0'}];
      %  PairList = [ PairList; {'img-manmade2-p-169t0','img-manmade5-p-282t0'}];
     %   PairList = [ PairList; {'img-manmade2-p-169t0','img-manmade5-p-282t0'}];
    %    PairList = [ PairList; {'img-manmade6-p-107t0','img-manmade7-p-313t0'}];
    %    PairList = [ PairList; {'img-manmade6-p-107t0','img-manmade7-p-313t0'}];
        PairList = [ PairList; {'img-manmade6-p-107t0','img-manmade7-p-343t0'}];
       % PairList = [ PairList; {'img-manmade7-p-313t0','img-manmade7-p-343t0'}];
FlagSetup;
%Flag.NewInput = 0;
%Flag.FlagDisp = 0;
%Flag.PoseMatches = 1;
%Flag.ReInference = 0;
%Flag.PoseMatchStor = 1;
%Flag.FindOcclu = 1;
%Flag.NewModelList = 1; % if 0 used the prestoraged tri-angulation data
%Flag.FlagFirstPair = Flag.NewModelList;

tic;
main( Fdir, PairList, WlrName, 1, Flag);
toc;

cd Batch/;

% Exit MATLAB
exit
EOF

# Display the output
cat Model.out
