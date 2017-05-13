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
function [SwissCoord]=LongLat2SwissCord(defaultPara, Geo)

% This function follow the instruction in u-blox document
% - GPS_Basics, and also wiki en.wikipedia.org/wili/Swiss_coordinate_system

Phi = ((Geo.Latitude*3600) - 169028.66)/10000
Gama = ((Geo.Longitude*3600) - 26782.5)/10000

SwissCoord(1) = 600072.37+(211455.93*Gama)-(10938.51*Gama*Phi)-(0.36*Gama*(Phi^2))-(44.54*(Gama^3));
SwissCoord(2) = 200147.07+(308807.95*Phi)-(3745.25*(Gama^2))-(76.63*(Phi^2))-(194.56*Phi*(Gama^2))+(119.79*(Phi^3));
SwissCoord(3) = (Geo.altitude - 49.55)+(2.73*Gama)+(6.94*Phi);

return;
