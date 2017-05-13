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
function [X] = Polar2Cartesian(defaultPara, Long, Lat, Alt, FlagConv)

% This function convert the Geographic coordinate system in polar coordinate to
% Cartesian coordinate system which set z in the north south pole and
% x pointing out from Africa see: http://en.wikipedia.org/wiki/Latitude_and_Longitude
% for more information
% implement using equation in
% http://www.colorado.edu/geography/gcraft/notes/datum/gif/llhxyz.gif
%
% Input:
%       defaultPara - ellip_equatorial_radius, defaultPara.ellip_polar_radius.
%       Long - Longitude in Decimal Degree (can be column vector)
%       Lat - Latitude in Decimal Degree
%       Alt - Altitude in Decimal Degree (above means see level)
%
% Return:
%       X - (x y z ) in Cartesian coordinate

if nargin <5
    FlagConv = false;
end
% convert degree to Radian
if FlagConv
    Long = Long/180*pi;
    Lat = Lat/180*pi;
end    

% calculate intermediant parameter
f_geo = (defaultPara.ellip_equatorial_radius - defaultPara.ellip_polar_radius)...
            /defaultPara.ellip_equatorial_radius; % flattening
tempf = (defaultPara.ellip_equatorial_radius^2 - defaultPara.ellip_polar_radius^2)...
            /(defaultPara.ellip_equatorial_radius)^2;
Eccen_sqr = 2*f_geo-f_geo^2; % eccentricity_squares
tt = sqrt(1 - Eccen_sqr*(sin(Lat)^2) );
NewN = defaultPara.ellip_equatorial_radius/sqrt( 1-(tempf*sin(Lat))^2);
N = defaultPara.ellip_equatorial_radius/tt;% radius of curvature in prime vertical

% main port of Cartesian2Polar
% x = (N+Alt).*cos(Lat).*cos(Long);
% y = (N+Alt).*cos(Lat).*sin(Long);
% z = ((N)*(1-Eccen_sqr)+Alt).*sin(Lat);%z = ((N)*(1)+Alt).*sin(Lat);
x = (NewN+Alt).*cos(Lat).*cos(Long);
y = (NewN+Alt).*cos(Lat).*sin(Long);
z = ((NewN)*(1-tempf)+Alt).*sin(Lat);


X = [x; y; z];
return;
