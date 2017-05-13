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
matlab -nojvm -nodisplay > ChainListModel"$Task"Corr.out << EOF

% Run your MATLAB commands inlin
cd ..
image3dstichingSetPath
Fdir ='/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/$Task';

% process the List
List = [];
% read in file from $ChainList
fp = fopen(['./Batch/' '$ChainList'],'r');
tline = fgetl(fp);	
while tline ~= -1
	List = [List; { tline}];
	tline = fgetl(fp);	
end
%        List = [ List; {  'IMG_0892', 'IMG_0895' }];
%        List = [ List; {  'IMG_0892', 'IMG_0896' }];
%        List = [ List; {  'IMG_0866', 'IMG_0867' }];

% run each List
for i = 1:size(List,1);

	ListTemp = List{i};	
	id = strfind( ListTemp,' ');
	id = [0 id length(ListTemp)+1];
	ChainList = [];
	Wrlname = '$Task';
	for i = 1:( length(id)-2)
		j = i +1;
		ChainList = [ ChainList; { ListTemp( ( id(i)+1 ):( id(i+1)-1)), ListTemp( ( id(j)+1 ):( id(j+1)-1))}];
		Wrlname = [ Wrlname '_' ListTemp( ( id(i)+1 ):( id(i+1)-1))];
	end
	Wrlname = [ Wrlname '_' ListTemp( ( id(j)+1 ):( id(j+1)-1))];

	%FlagSetup;

	FlagSetupTestPairNew;
	% ======= Flag Changes
	Flag.FlagImgInfoLoadPreStorage = 3;
	Flag.PoseMatches = 1;
	Flag.ReInference = 1;
	Flag.Refinement = 1;

	% Refinement
	Flag.FlagPreloadOccluDetectMatches = 1;
	Flag.LoadStoragedSurfOccluMatches = 1;
	Flag.FlagPreloadOccluDetectRay = 1;
	Flag.loadStoragedPairNew = 1;
	Flag.LoadCorrMatches = 0;
	% ====================


	Stiching3dParameterSetup;
	% ======= Parameter Changes ======
	defaultPara.AbsThre = $AbsThre;
	defaultPara.RatioThre = $RatioThre;
	% ================================
	tic;
	main( Fdir, ChainList, Wrlname, defaultPara);
	toc;
end

cd Batch/;

% Exit MATLAB
exit
EOF

# Display the output
cat ChainListModel"$Task"Corr.out

