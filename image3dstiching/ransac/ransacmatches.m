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
function [F, inliers, fail] = ransacmatches(defaultPara, f1, f2, matches, I1, I2, disp)
% Computes the fundamental matrix and inlier matches using ransac.  Points
% are sampled non-uniformly in order to prefer more matches that are spread
% across the image.  Otherwise the algorithm is standard.
% input:  f1          - x,y coordinates of all feature frames in image 1
%         f2          - same for image 2
%         matches     - 2 by nummatches array specifying the initial set of
%                       possible matches between f1 and f2
%         I1/I2       - optional images to display
%         disp        - if true, display the matches found when done.

x1 = [];
x2 = [];

nmatches = size(matches, 2)
for i=1:nmatches
    x1(i, 1:2) = f1(1:2, matches(1, i));
    x2(i, 1:2) = f2(1:2, matches(2, i));
end

% calculate distances
% d1 is the sum of the squares of the distance from each point
% in im1 to every other point
% will use dist (the normalized avg distance for that match)
% to weight the sampling algorithm for ransac
d1=[];
d2=[];
for i=1:nmatches
    d1(i) = sum(sum( ((ones(nmatches, 1) * x1(i, 1:2)) - x1).^2));
    d2(i) = sum(sum( ((ones(nmatches, 1) * x2(i, 1:2)) - x2).^2));
end
dist=(d1+d2)/sum(sum(d1+d2));

% Assemble homogeneous feature coordinates for fitting of the
% fundamental matrix, note that [x,y] corresponds to [col, row]
x1 = [x1'; ones(1, length(x1))]; %[m1(2,:); m1(1,:); ones(1,length(m1))];
x2 = [x2'; ones(1, length(x1))]; %[m2(2,:); m2(1,:); ones(1,length(m1))];    
    
t = .002;  % Distance threshold for deciding outliers
t_more = .001;  % Distance threshold for deciding outliers (MS added)

[F, inliers, fail] = ransacfitfundmatrix(defaultPara, x1, x2, t, zeros(2), zeros(2), dist, 1, 1, 0);
% [inliers_more] = ransacPrune(F, x1, x2, t_more, 1, dist); % MS added

if disp
    figure(2) ; clf ;
    plotmatches(I1,I2,f1, f2,matches(:, inliers), 'Stacking', 'v'); %, 'Interactive', 1) ;
    vgg_gui_F(I1, I2, F');
end
