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
function [Info] = ReadLookAtInfo(out, Info, Index)

% Read info of longitude latitude altitude tilt heading roll
ptr1 = findstr(out, '<');
if numel(ptr1) == 0
   disp('check'); 
end    
ptr2 = findstr(out, '>');
Target = out( ( ptr1(1)+1): ( ptr2(1)-1));
StartPtr = ptr2(1) + 1;
switch lower(Target)
	case 'longitude'
  		EndPtr = findstr(out, '</longitude>') - 1;
		Info(Index).Geo.Longitude = str2num(out(StartPtr:EndPtr));
	case 'latitude'
  		EndPtr = findstr(out, '</latitude>') - 1;
		Info(Index).Geo.Latitude = str2num(out(StartPtr:EndPtr));
	case 'altitude'
  		EndPtr = findstr(out, '</altitude>') - 1;
		Info(Index).Geo.altitude = str2num(out(StartPtr:EndPtr));
	case 'tile'
  		EndPtr = findstr(out, '</tilt>') - 1;
		Info(Index).Rw(1) = str2num(out(StartPtr:EndPtr));
	case 'heading'
  		EndPtr = findstr(out, '</heading>') - 1;
		Info(Index).Rw(2) = str2num(out(StartPtr:EndPtr));
	case 'roll'
  		EndPtr = findstr(out, '</roll>') - 1;
		Info(Index).Rw(3) = str2num(out(StartPtr:EndPtr));
end
