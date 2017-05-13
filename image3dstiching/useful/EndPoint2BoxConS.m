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
function [ Rc1_2 ConS1_2 RoughConS1_2 tempConS1_2] = EndPoint2BoxConS(defaultPara, H, V, x1_2Max, x1_2Min, ExpandAlongEpipoLineFlag)

% [ Rc1_2 ConS1_2 RoughConS1_2] = EndPoint2BoxConS(defaultPara, H, V, x1_2Max, x1_2Min)
% data structure
% Rc1_2 - 4 by Num of x1_2Max = size(x1_2Max,2)

Ratio2RegularMatches = 0.1;
ExpandRatio = 1.25;

% if ExpandAlongEpipoLineFlag == 1
% expand the search space along the epipolar line
if ExpandAlongEpipoLineFlag
   UnitVector = x1_2Max - x1_2Min;
   Len = norm(UnitVector);
   UnitVector = UnitVector/Len;
   Center = (x1_2Max + x1_2Min)/2;
   % Calculate new x1_2Max and x1_2Min
   x1_2Max = Center + UnitVector*(Len/2)*ExpandRatio;
   x1_2Min = Center - UnitVector*(Len/2)*ExpandRatio;
end

% used the same concept in ../match/CalMatchSearchRegin.m

	tRAN = x1_2Max - x1_2Min;    	
	
	Ptr = tRAN(1,:) < 0;
	%theta_z = -atan(tRAN(2,:)./tRAN(1,:));
	theta_z = atan(tRAN(2,:)./tRAN(1,:));
	theta_z(Ptr) = theta_z(Ptr)+pi;
	Rc1_2 = [cos(-theta_z); -sin(-theta_z); sin(-theta_z); cos(-theta_z)];
	ConS1_2(1:2,:) = x1_2Min;
	ConS1_2(3,:) = sum( Rc1_2(1:2,:).*tRAN,1);% also may check sum( Rc1_2(3:4,:).*tRAN,1) close up to zeros to verify correct or not
%     Re1_2 = sum( Rc1_2(3:4,:).*tRAN,1)
	
	ConS1_2(4,:) = defaultPara.VertVar*max(H,V)*Ratio2RegularMatches;
        % [x_origon; y_origin; x_bound, y_bound(defaultPara.VertVar*max(H,V))]

        % generate plot point
	tempConS1_2 = [ ConS1_2(3,:); ConS1_2(4,:); ...
                        zeros(1,size(ConS1_2,2)); ConS1_2(4,:);...
			zeros(1,size(ConS1_2,2)); -ConS1_2(4,:);...
			ConS1_2(3,:); -ConS1_2(4,:)];
	%Rc1_2_transpose = [cos(theta_z); sin(theta_z); -sin(theta_z); cos(theta_z)];
	Rc1_2_transpose = [cos(theta_z); -sin(theta_z); sin(theta_z); cos(theta_z)];
	tempConS1_2(1:2,:) = [sum( Rc1_2_transpose(1:2,:).*tempConS1_2(1:2,:), 1); ...
			      sum( Rc1_2_transpose(3:4,:).*tempConS1_2(1:2,:), 1)];	
	tempConS1_2(3:4,:) = [sum( Rc1_2_transpose(1:2,:).*tempConS1_2(3:4,:), 1); ...
			      sum( Rc1_2_transpose(3:4,:).*tempConS1_2(3:4,:), 1)];
	tempConS1_2(5:6,:) = [sum( Rc1_2_transpose(1:2,:).*tempConS1_2(5:6,:), 1); ...
			      sum( Rc1_2_transpose(3:4,:).*tempConS1_2(5:6,:), 1)];
	tempConS1_2(7:8,:) = [sum( Rc1_2_transpose(1:2,:).*tempConS1_2(7:8,:), 1); ...
			      sum( Rc1_2_transpose(3:4,:).*tempConS1_2(7:8,:), 1)];
	tempConS1_2 = tempConS1_2 + repmat( x1_2Min, 4, 1);
    RoughConS1_2 = [ min( tempConS1_2([ 1 3 5 7],:), [], 1); max( tempConS1_2([ 1 3 5 7],:), [], 1);...
                    min( tempConS1_2([ 2 4 6 8],:), [], 1); max( tempConS1_2([ 2 4 6 8],:), [], 1)];

	% check if tRAN is zero
	MaskZeros = tRAN(1,:) == 0;
	Rc1_2(:, MaskZeros) = 0;
	ConS1_2(:,MaskZeros) = 0;
	RoughConS1_2(:,MaskZeros) =  -1;
return;
