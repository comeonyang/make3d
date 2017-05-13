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
function [ P, X, Outliner] = ProjectionFactorization_missing_data...
		( ARes, BRes, x_im, inc, depth)
% This funciton given the match in image coordinate and the initial depth
% and estimate the Camera matrix and the 3d position of the match points

% setup parameters
Min_prjective_depth_change_ratio = 1e-6;
ProjErrorRatio = 0.0005;
No_Camera = size(depth,3); 
No_Points = size(x_im,2)
MaxCount = 100;
TotalOutlinerPercentage = 50;
IterOutlinerPercentage = TotalOutlinerPercentage/MaxCount;

% building matrix 
lamda = [];
Xim = [];
INC = [];
for i = 1:No_Camera
    lamda = [ lamda; repmat(depth(:,:,i), 3, 1)]; 
    Xim_miss = [ Xim; [x_im(:,:,i); ones( 1, No_Points )]];
    INC = [INC; repmat(inc(:,:,i),3,1)];
end

count = 1
Outliner = [];
while count < MaxCount
    
    % normalized the depth and Xim
    RowNormalizeFactor = sum(lamda,2);
    RowNormalizeCount = sum(lamda ~=0, 2);
    RowNormalizeFactor = RowNormalizeFactor./RowNormalizeCount;
    RowNormalizeFactor = sqrt(1./RowNormalizeFactor);
    %lamdaNormal = lamda.*repmat( RowNormalizeFactor, 1, No_Points);
    ColumnNormalizeFactor = sum(lamda,1);
    ColNormalizeCount = sum(lamda ~=0, 1);
    ColNormalizeFactor = ColNormalizeFactor./ColNormalizeCount;
    ColNormalizeFactor = sqrt(1./ColNormalizeFactor);
    %lamdaNormal = lamdaNormal.*repmat( ColumnNormalizeFactor,No_Camera*3,1);
    lamdaNormal = RowNormalizeFactor*ColNormalizeFactor;

    % First Complete missing data
    [e,Xim,s] = rankrsfm(lamdaNormal.*Xim_miss,INC);

    % use estimated depth for missing data
    for i = 1:No_Camera
 	lamda((3*i-2):(3*i),find(~inc(:,:,i))) = ...
		repmat( Xim(3*i,find(~inc(:,:,i)))./lamdaNormal(3*i,find(~inc(:,:,i))), 3, 1);	
    end

    % normalized the depth and Xim again with new lamda
    RowNormalizeFactor = 1./( sqrt( sum(lamda.^2,2)));
    lamdaNormal = lamda.*repmat( RowNormalizeFactor, 1, No_Points);
    ColumnNormalizeFactor = sqrt(3)./( sqrt( sum(lamdaNormal.^2,1)));
    lamdaNormal = lamdaNormal.*repmat( ColumnNormalizeFactor,No_Camera*3,1);

    % Solving SVD
    [ U S V] = svds(lamdaNormal.*Xim,4);
    P = U*S;
    X = V';

    % new lamda
    M = P*X;
%    M = M./repmat( RowNormalizeFactor, 1, No_Points);
%    M = M./repmat( ColumnNormalizeFactor,No_Camera*3,1);
    M = M./lamdaNormal;
    ThirdSample = (1:No_Camera)*3;
    TempLamda = M( ThirdSample, :);
 
    % calculate the geometric differency of measured and estimated projection points on the image plane
    NewLamda = TempLamda( reshape( repmat( 1:No_Camera, 3, 1), [], 1),:);
    x_im_est = M./NewLamda;
    xIm = Xim( reshape( [ ThirdSample-2; ThirdSample-1], [], 1), :);
    xImEst = x_im_est( reshape( [ ThirdSample-2; ThirdSample-1], [], 1), :);
    DiffProj =  xIm - xImEst;
    EclideanError = norms( reshape(DiffProj, 2, []));
    InlinerUpperBound = prctile( EclideanError, 100 - IterOutlinerPercentage);
    OutlinerPtr = EclideanError > InlinerUpperBound;
    OutlinerPtr = reshape( OutlinerPtr, 2, []);
    OutlinerPtr = sum( OutlinerPtr,1);
    figure(15);
    hist( EclideanError,100);
    AAvgDiffProj(count) = mean(norms(DiffProj(1:2,:)));
    BAvgDiffProj(count) = mean(norms(DiffProj(3:4,:)));
    if AAvgDiffProj(count) < norm( ARes*ProjErrorRatio) && BAvgDiffProj(count) < norm( BRes*ProjErrorRatio)
       break;
    end    
    
    % Check stop
    Ratio(count) = max( max( abs((NewLamda( ThirdSample, :) - lamda( ThirdSample, :))./lamda( ThirdSample, :))));
    if Ratio(count) < Min_prjective_depth_change_ratio
       break;
    end
    count = count +1;
    lamda = NewLamda;
    
    % =============remove Outliner =============
    Outliner = [ Outliner find(OutlinerPtr~=0)];
    lamda(:, OutlinerPtr~=0) = [];
    Xim(:, OutlinerPtr~=0) = [];
    No_Points = size( Xim,2);
    % ==========================================
end
No_Points
count
AAvgDiffProj(count-1)
BAvgDiffProj(count-1)
figure(100); hist( DiffProj(:), 1000);
figure(10); scatter(xImEst(1,:),xImEst(2,:),4,ones(1,size(xImEst,2)));
hold on; scatter(xIm(1,:),xIm(2,:),4,0.5*ones(1,size(xIm,2)));hold off;
figure(11); scatter(xImEst(3,:),xImEst(4,:),4,ones(1,size(xImEst,2)));
hold on; scatter(xIm(3,:),xIm(4,:),4,0.5*ones(1,size(xIm,2)));hold off;

% restorge P X
 P = P./repmat( RowNormalizeFactor, 1, 4);
 X = X./repmat( ColumnNormalizeFactor, 4, 1);
%  X = X./repmat(X(4,:), 4,1);

return;
