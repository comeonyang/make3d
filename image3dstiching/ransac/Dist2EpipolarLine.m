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
function [d] = Dist2EpipolarLine(F, x)

% This function calculate Both dist of x1 to l1 and x2 to l2
% Which are both the distance to the epipolar line

% Initialize variables
Nmatches = size(x,2);
if size(x,1) == 4
	x1 = [ x(1:2,:); ones(1, Nmatches)];    % Extract x1 and x2 from x
	x2 = [ x(3:4,:); ones(1, Nmatches)];
else
	x1 = x(1:3,:);    % Extract x1 and x2 from x
	x2 = x(4:6,:);
end

%
x2tFx1 = zeros(1, Nmatches);
for n = 1:Nmatches
	x2tFx1(n) = x2(:,n)'*F*x1(:,n);
end

Fx1 = F*x1;
Ftx2 = F'*x2;

% evaluate distance in pixels
d(1,:) =  x2tFx1 ./ ...
	sqrt(Fx1(1,:).^2 + Fx1(2,:).^2 );
d(2,:) =  x2tFx1 ./ ...
	sqrt(Ftx2(1,:).^2 + Ftx2(2,:).^2);
d(3,:) =  x2tFx1.^2 ./ ...
      (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
return; 
