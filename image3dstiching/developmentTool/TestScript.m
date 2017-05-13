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
% variable
Task = 'QuadCorner';
PairListFlag = 1;
PairList = 'PairList.txt';
AbsThre = 0.2;
RatioThre = 0.6;

% image3dstichingSetPath
Fdir =['/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/' Task];

% process the List
if PairListFlag
        List = [];
        % read in file from $PairList
        fp = fopen(['./Batch/' PairList],'r');
        tline = 1;      
        tline = fgetl(fp);
        while tline ~= -1                
                id = strfind(tline,' ');
                List = [ List; { tline(1:( id-1 )), tline(id:end)} ];
                tline = fgetl(fp);
        end
else
        file = dir([ Fdir '/jpg/*.jpg']);
        List = [];
        for i = 1:length(file)
                for j = 1:length(file)
                        if i<j
                                List = [List; { strrep(file(i).name,'.jpg',''), strrep(file(j).name,'.jpg','')}];
                        end
                end
        end
end

% run each List
for i = 1:size(List,1);
        WlrName = [ Task '_' List{i,1} '_' List{i,2}];
        PairList = [];
        PairList = [ PairList; {List{i,1},List{i,2}}];

        %FlagSetup;

        FlagSetupTestPairNew;
        % ======= Flag Changes
        Flag.FlagImgInfoLoadPreStorage = 4;
        Flag.PoseMatches = 1;
        Flag.ReInference = 1;
        Flag.Refinement = 0;

        % Refinement
        Flag.FlagPreloadOccluDetectMatches = 0;
        Flag.LoadStoragedSurfOccluMatches = 0;
        Flag.FlagPreloadOccluDetectRay = 0;
        Flag.LoadCorrMatches = 0;
        % ====================


        Stiching3dParameterSetup;
        % ======= Parameter Changes ======
        defaultPara.AbsThre = AbsThre;
        defaultPara.RatioThre = RatioThre;
        % ================================
        tic;
        main( Fdir, PairList, WlrName, defaultPara);
        toc;
end
