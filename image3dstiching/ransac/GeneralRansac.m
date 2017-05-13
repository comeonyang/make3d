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
function [F, inliers, NewDist, fail, ind]=GeneralRansac(defaultPara, f1, f2, matches, D1, D2, PriorDist, s);

        nmatches = size(matches, 2);

	if nargin < 7
		PriorDist = ones(1,nmatches);
                s = 8;
	elseif nargin <8
		s = 8;
	end

    for i=1:nmatches
            x1(1:2, i) = f1(1:2, matches(1, i));
            x2(1:2, i) = f2(1:2, matches(2, i));
    end
    
    %t = .0002;  % Distance threshold for deciding outliers
    t = .00002;  % Distance threshold for deciding outliers % July 30 Min added worked better

    % Initialize distribution to uniform
    dist = ones(nmatches,1); % using uniform dist gives better result
    dist = dist./sum(dist);
	
    % Condition on Prior Dist
	dist = dist.*PriorDist; % NOTICE!! The new dist is not normalized


    % First Step Ransac to define dist from learned Depth
    [F, inliers, NewDist, fail, ind] = ransacfitfundmatrix(defaultPara, x1, x2, t, D1, D2, dist, 0, 1, 0, s);

return;
