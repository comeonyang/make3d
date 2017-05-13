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
function [F, newinliers, fail] = RansacNonlinearReFinement(f1, f2, matches, F, inliers, I1, I2, disp)
% first re-estimates F based on a non-linear method, then projects each
% point into 3D space and removes matches that are either very close or
% very far away from camera 1

numMatches = size(matches, 2);
numInliers = size(inliers, 2);
percentInliers = numInliers/numMatches;

fail = 0;
if (percentInliers < 0.2) | (numInliers < 20)
    fail = 1;
    newinliers = inliers;
else
    x1 = [f1(:, matches(1, inliers))];
    x2 = [f2(:, matches(2, inliers))];
    x = [x1' x2'];
    m3=1; %256

    % re-estimate the fundamental matrix using a non-linear method
    [f, f_sq_errors, n_inliers,inlier_index,F] = torr_estimateF( x, m3, [], 'non_linear', 0, F);

    % find the 3D positions of the points
    X = ThreeDPoints(F, size(I1), size(I2), ...
        f1(:, matches(1, inliers)), f2(:, matches(2, inliers)));%F'??

    size(X)

    % find the distance to camera1, sort them, find the size of gaps
    % between points
    dist = sqrt(sum(X.^2));
    sdist = sort(dist);
    spans = sdist - [0 sdist(1:end-1)];

    % find the start/end dist of the middle 50% of points
    numpts = length(sdist)
    startLind = ceil(numpts*0.25);
    startHind = ceil(numpts*0.75);
    startL = sdist(startLind);
    startH = sdist(startHind);

    % don't allow spans greater than 20 times the average span of the
    % middle 50% of points
    maxSpan = 20*(startH - startL)/numpts;
    badspansH = find(spans>maxSpan);
    badspansL = find(spans>maxSpan/2);
    [lowidx, temp] = max(badspansL(badspansL<startLind));
    if isempty(lowidx)
        lowidx=1;
    end
    [highidx, temp] = min(badspansH(badspansH>startHind));
    if isempty(highidx)
        highidx=length(sdist);
    end
    lowbound = sdist(lowidx)
    highbound = sdist(highidx)
    keep = (dist>lowbound).*(dist<highbound);
    idx = find(keep);
    rejectIdx=find(1-keep);

    % get the new inliers
    newinliers = inliers(idx);
    numNewInliers = length(newinliers)
    newoutliers = inliers(rejectIdx);

    %get new fund matrix
    %x1 = [f1(:, matches(1, newinliers))];
    %x2 = [f2(:, matches(2, newinliers))];
    %x = [x1; ones(1, length(x1)); x2; ones(1, length(x1))];
    %F = fundmatrix(x);

    %re-evaluate all points
    %x1 = [f1(:, matches(1, :))];
    %x2 = [f2(:, matches(2, :))];
    %x = [x1; ones(1, length(x1)); x2; ones(1, length(x1))];
    %refitinliers = getInliers(x, F, 0.00005);

    if (disp)
%         figure(3)
%         hist(sdist, 100)
       figure;
       vgg_gui_F(I1, I2, F');

        figure;
        clf;
        plotmatches(I1,I2,f1, f2,matches(:, newinliers), 'Stacking', 'v', 'Interactive', 2) ;
%         figure(3)
%         clf;
%         plotmatches(I1,I2,f1, f2,matches(:, inliers), ...
%             'Stacking', 'v', 'Interactive', 1, 'Dist', dist) ;
        %plotmatches(I1,I2,f1, f2,matches(:, refitinliers), 'Stacking', 'v', 'Interactive', 1) ;
    end
end

