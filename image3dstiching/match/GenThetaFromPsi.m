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
function [theta]= GenThetaFromPsi(x1, x2, psi);

% This function is called by Matches2ViewpointTransform
% give matches point and psi samples
% calcuated each theta for all matches and all samples

NumSample = length(psi);
NumMatches = length(x1);

% first step rotation * x1
R11 = permute( cos(psi), [ 1 3 2]);
R13 = permute( sin(psi), [ 1 3 2]);
Newx1 =[ repmat( x1(1,:), [1 1 NumSample]).*repmat( R11, [1 NumMatches 1]) + repmat( R13, [1 NumMatches 1]);...
	 repmat( x1(2,:), [1 1 NumSample]);...
         repmat( x1(1,:), [1 1 NumSample]).*repmat( -R13, [1 NumMatches 1]) + repmat( R11, [1 NumMatches 1])]; %(3 x NumMatches x NumSample)

% Last step find theta for each entries of Newx1
Interm = (Newx1(1,:,:).*repmat( x2(2,:), [1 1 NumSample]) - Newx1(2,:,:).*repmat( x2(1,:), [1 1 NumSample]))./...
	 (Newx1(3,:,:).*repmat( x2(2,:), [1 1 NumSample]) + Newx1(2,:,:));

theta = atan( permute( Interm, [ 2 3 1]));

