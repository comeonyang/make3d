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
function [] = NewVrmlByPick(Fdir, Wrlname, ImgList, OutDir)

% function generate new inline VRML uing Wrlname nad ImgList in Outdir folder

% method keep deleting the .wrl not in the ImgList
NewModel = cell2str(ImgList);
SearchList = ImgList;
Path = [Fdir '/' Wrlname '.wrl'];
fp = fopen( [ OutDir '/' NewModel '.wrl'], 'w');
fprintf(fp, '#VRML V2.0 utf8\n');
fclose(fp);
while length(SearchList)~=0
	fp = fopen( Path, 'r');
    line = 1;
    out = 1;
	while out ~= -1
        out = fgets(fp);
		if ~isempty( findstr(out,SearchList{1}))
%			[ 'sed -n ' num2str(line-6) ',' num2str(line+3) 'p ' '< ' Path ' >> ' OutDir '/' NewModel '.wrl']
			system([ 'sed -n ' num2str(line-6) ',' num2str(line+3) 'p ' '< ' Path ' >> ' OutDir '/' NewModel '.wrl']);
			break;
        end
        line = line +1;
	end
	SearchList(1) = [];
    fclose(fp);
end

return;
