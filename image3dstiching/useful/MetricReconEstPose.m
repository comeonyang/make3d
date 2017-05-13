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
function [R_ground_constrain T_mr D1_modified D2_modified] = MetricReconEstPose(P, R, T, T1, T2, D1, D2, x1, x2, defaultPara);

% This function reconstruct the R_mr and T_mr from the P( camera matrix)
% which most close to the original estimated R and T
    GroundStaticWeight = 50;
    ops = sdpsettings('solver','sedumi','verbose',1);
    
    P = [inv(T1)*P(1:3,:);...
         inv(T2)*P(4:6,:)];
    
    P1 = P(1:3,:);
    P2 = P(4:6,:);
    [U S V] = svd(P1);
    H_last = V(:,end);
    H_4by3 = P1\defaultPara.InrinsicK1;
    
    K2simple_square = sdpvar(3,1);
    F = set(diag(K2simple_square) >=0);
    sol = solvesdp(F , norm( P2*H_4by3*H_4by3'*P2' - diag(K2simple_square),'fro'), ops);
    K2simple_square = double(K2simple_square);
    K2simple = sqrt(K2simple_square);
    R_est = diag( 1./K2simple)*P2*H_4by3;
%     R_est = P2*H_4by3;
    T_est = P2*H_last;
    
    % Reconstruct Rotation matrix by penalize rotation in other than z axis
    R_ground_constrain = sdpvar(3,3);
    sol =solvesdp( [], norm( R_ground_constrain - R_est,'fro') + GroundStaticWeight*norm([0 1 0]' - R_ground_constrain*[0 1 0]'), ops);
    R_ground_constrain = double( R_ground_constrain);
    R_ground_constrain =  R_ground_constrain * (R_ground_constrain'*R_ground_constrain)^(-.5);

    % Reconstruct the translation matrix by using Mono-depth infomation
    Ray1 = inv(defaultPara.InrinsicK1)*x1;
    Ray2 = inv(defaultPara.InrinsicK2)*x2;
    T = sdpvar(3,1);
    D1_modified = sdpvar(1,length(D1));
    D2_modified = sdpvar(1,length(D2));
    D2scale = sdpvar(1);
    sol =solvesdp( [], norm( reshape( (Ray1.*repmat( D1_modified, 3,1))' - ...
                    ( R_ground_constrain'*( Ray2.*repmat( D2_modified, 3, 1) ) + repmat(T,1,length(D1)) )' ,1,[]),1)+...
                    norm(D1_modified - D1,1) + norm(D2_modified - D2scale*D2,1), ops);
    T_mr = double(T);
    D1_modified = double(D1_modified);
    D2_modified = double(D2_modified);
    D2scale = double(D2scale);

   return; 
