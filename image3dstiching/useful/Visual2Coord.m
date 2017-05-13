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
function []=Visual2Coord(R,T);

% This function plot the relation of the two coordinate

% initial parameter
Axis_length = 10;
z_length = 50;

% original
center =zeros(3,1);
top = Axis_length*[0 1 0]';
left = Axis_length*[1 0 0]';
front = z_length*[0 0 1]';

% transform
Tcenter = center + T;
Ttop = R*top + T;
Tleft = R*left + T;
Tfront = R*front + T;

figure;
scatter3([center(1) Tcenter(1)],[center(3) Tcenter(3)],[center(2) Tcenter(2)]);
hold on;
line([center(1) left(1)]', [center(3) left(3)]', [center(2) left(2)]', 'Color' ,'g');
line([Tcenter(1) Tleft(1)]', [Tcenter(3) Tleft(3)]', [Tcenter(2) Tleft(2)]', 'Color' ,'g');
line([center(1) top(1)]', [center(3) top(3)]', [center(2) top(2)]', 'Color' ,'r');
line([Tcenter(1) Ttop(1)]', [Tcenter(3) Ttop(3)]', [Tcenter(2) Ttop(2)]', 'Color' ,'r');
line([center(1) front(1)]', [center(3) front(3)]', [center(2) front(2)]', 'Color' ,'y');
line([Tcenter(1) Tfront(1)]', [Tcenter(3) Tfront(3)]', [Tcenter(2) Tfront(2)]', 'Color' ,'y');
axis equal;
