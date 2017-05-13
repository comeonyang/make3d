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
function [ImgInfo]=ExtractGoogleEarthInfo(ImgInfo, FPath)

% This function extract the info in Kml file
% - June 26:
%		It extract longitude and latitude and heading info

fp = fopen(FPath,'r');
out = fgets(fp);
while out ~=-1
	if ~isempty( findstr(out, '<Placemark>'))
		out = fgets(fp);
		ptr = findstr(out, '<name>');
		if ~isempty(ptr)
			StartPtr = ptr(1) + 6;
			EndPtr = findstr(out, '</name>') - 1;	
			name = out(StartPtr:EndPtr);
			ImgIndex = ImgInfoIndexFromName(ImgInfo, name);
		
			% Start extracting info in <LootAt> and </LookAt>
			out = fgets(fp);
			ptr1 = findstr(out, '</LookAt>');
			ptr2 = findstr(out, '</Placemark>');
			while isempty(ptr1) && isempty(ptr2)
				ImgInfo = ReadLookAtInfo(out, ImgInfo, ImgIndex);
                out = fgets(fp);
                ptr1 = findstr(out, '</LookAt>');
                ptr2 = findstr(out, '</Placemark>');
			end
		end
	end
	out = fgets(fp);
end
