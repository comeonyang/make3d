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
function [ F, inliers, fail]=RansacOnPairMatches(defaultPara, f1, f2, matches, I1, I2, Depth1, Depth2, disp)

% This function Use the Ransac algorithm to generate models 
% with Fundamental matrix and inlier matches
% that suit the matches most
% Input:
%          defaultPara - useful default parameters (like, camera intrinsic
%          matrix)
%         f1          - x,y coordinates of all feature frames in image 1
%         f2          - same for image 2
%         matches     - 2 by nummatches array specifying the initial set of
%                       possible matches between f1 and f2
%         I1/I2       - optional images to display
%         Depth1/2    - depth imformation to support more accurate ransac
%         disp        - if true, display the matches found when done.
% See Also:
%         ransacfitfundmatrix.m (in kovesi)

x1 = [];
x2 = [];

nmatches = size(matches, 2);
for i=1:nmatches
    x1(i, 1:2) = f1(1:2, matches(1, i));
    x2(i, 1:2) = f2(1:2, matches(2, i));
end

% Assemble homogeneous feature coordinates for fitting of the
% fundamental matrix, note that [x,y] corresponds to [col, row]
x1 = [x1'; ones(1, length(x1))]; %[m1(2,:); m1(1,:); ones(1,length(m1))];
x2 = [x2'; ones(1, length(x1))]; %[m2(2,:); m2(1,:); ones(1,length(m1))];    
X = [inv(defaultPara.InrinsicK1)*x1;...
     inv(defaultPara.InrinsicK2)*x2 ];
t = .002;  % Distance threshold for deciding outliers

% Initialize distribution to uniform
dist = ones(nmatches,1); % using uniform dist gives better result
dist = dist./sum(dist);

% First Step Ransac to define dist from learned Depth
[F, inliers, NewDist, fail] = ransacfitfundmatrix(defaultPara, x1, x2, t, Depth1, Depth2, dist, 1, disp, 0);
if disp
%    figure(2) ; clf ;
%    plotmatches(I1,I2,f1, f2,matches(:, inliers), 'Stacking', 'v', 'Interactive', 0) ;
    figure(2) ; clf ;
    plotmatches(I1,I2,x1(1:2,:), x2(1:2,:), [inliers; inliers], 'Stacking', 'v', 'Interactive', 3) ;
end 

% Estimated Distribution from the F and Depth info
% assume the know the camera intrinsic parameter of both camera
[x1, x2, Depth1, Depth2, X, lamda1, lamda2] = RemoveOutlier(x1, x2, Depth1, Depth2, X, [], [], inliers);
[R, T, lamda1, lamda2, inlierPosD] = EstPose( ...
       defaultPara.InrinsicK2'*F*defaultPara.InrinsicK1, X);
if disp
    figure(3) ; clf ;
    plotmatches(I1,I2,x1(1:2,:), x2(1:2,:), repmat(setdiff(1:size(X,2), inlierPosD), 2, 1), 'Stacking', 'v', 'Interactive', 0) ;
end 
% Calculate EstDepMatchDist
[x1, x2, Depth1, Depth2, X, lamda1, lamda2] = RemoveOutlier(x1, x2, Depth1, Depth2, X, lamda1, lamda2, inlierPosD);
[dist, inlierThreDist] = EstDepMatchDist(X, R, T, Depth1, Depth2, lamda1', lamda2', disp);
if disp
    figure(4) ; clf ;
    plotmatches(I1,I2,x1(1:2,:), x2(1:2,:),  repmat(setdiff(1:size(X,2),inlierThreDist),2,1), 'Stacking', 'v', 'Interactive', 3) ;
end

% Second Ransac with new Distribution
[x1, x2, Depth1, Depth2, X, lamda1, lamda2] = RemoveOutlier(x1, x2, Depth1, Depth2, X, lamda1, lamda2, inlierThreDist);
[F, inliers, NewDist, fail] = ransacfitfundmatrix(defaultPara, x1, x2, t, Depth1, Depth2, dist, 1, disp, 0);

if disp
    figure(5) ; clf ;
    plotmatches(I1,I2,x1(1:2,:), x2(1:2,:),  repmat(inliers,2,1), 'Stacking', 'v', 'Interactive', 3) ;
%    vgg_gui_F(I1, I2, F');
end 
