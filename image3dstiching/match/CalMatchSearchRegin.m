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
function [ Rc1_2, Rc2_1, ConS1_2, ConS2_1, RoughConS1_2, RoughConS2_1]=CalMatchSearchRegin(defaultPara, R, T, I1, I2, x1, x2, D1, D2, FlagRotate, FlagDisp);

% This function Calculate the Constrain of SurfMatch search space give
% Estimated Rotation and Translation matrix  and depth
% in the reference camera coordinate

% initialize Parameter
NegativeDepthTolerence = defaultPara.NegativeDepthTolerence;
MaxRatio =defaultPara.MaxRatio;%300;
MinRatio = 1/MaxRatio;
R1_2 = R(1:3,:);
R2_1 = R(4:6,:);
T1_2 = T(1:3,:);
T2_1 = T(4:6,:);
K1 = size(x1,2);
K2 = size(x2,2);
x1 = [x1; ones(1,K1)];
x2 = [x2; ones(1,K2)];
%D1 = D(1,:);
%D2 = D(2,:);

% I1 project on I2 ==========================================% use zeros depth constraon to restrict the searching area more
[V H] = size(I2);
MaxD1 = min(D1*MaxRatio,defaultPara.FarestDist);
MinD1 = max(D1*MinRatio,defaultPara.Closestdist);
% calculate the prespective depth gives 0 prejective depth in camera 2
x1_2R_z = R1_2(3,:)*inv(defaultPara.InrinsicK1)*(x1); % it's possible that x1_2R_z is zero
mark =  x1_2R_z ~=0;
Depth_zero_ca2 = zeros(1, length(x1_2R_z));
Depth_one_ca2 = zeros(1, length(x1_2R_z));
Depth_zero_ca2(mark) = T1_2(3)./x1_2R_z(mark); % 0 (zero depth in camera 2 coordinate) = A*depth + T(3); solve depth
Depth_one_ca2(mark) = (1-T1_2(3))./x1_2R_z(mark); % 1 (zero depth in camera 2 coordinate) = A*depth + T(3); solve depth

% use zeros depth constraon to restrict the searching area more
mark0_1 = Depth_zero_ca2 >= Depth_one_ca2;
mark1_0 = ~mark0_1;

tempDepth_one_ca2 = Depth_one_ca2(mark0_1);
tempDepth_zero_ca2 = Depth_zero_ca2(mark0_1);
tempMax = MaxD1(mark0_1);
tempMin = MinD1(mark0_1);
markMaxBigerOneDepthCa2 = tempMax > tempDepth_one_ca2;
markMinBigerOneDepthCa2 = tempMin > tempDepth_one_ca2; % markMinBigerOneDepthCa2 : is one kind of outlier
tempMax(markMaxBigerOneDepthCa2) = tempDepth_one_ca2(markMaxBigerOneDepthCa2);
tempMax(markMinBigerOneDepthCa2) = tempDepth_zero_ca2(markMinBigerOneDepthCa2);%NaN; % indicator of outliers
tempMin(markMinBigerOneDepthCa2) = tempDepth_one_ca2(markMinBigerOneDepthCa2); % NaN indicator of outliers
MaxD1(mark0_1) = tempMax;
MinD1(mark0_1) = tempMin;
MaxD1(mark) = min(D1(mark)*MaxRatio,defaultPara.FarestDist);%Min529 new
MinD1(mark) = max(D1(mark)*MinRatio,defaultPara.Closestdist);%Min529 new


tempDepth_one_ca2 = Depth_one_ca2(mark1_0);
tempDepth_zero_ca2 = Depth_zero_ca2(mark1_0);
tempMax = MaxD1(mark1_0);
tempMin = MinD1(mark1_0);
markMaxBigerOneDepthCa2 = tempMax < tempDepth_one_ca2; % markMaxBigerOneDepthCa2 : is one kind of outlier
markMinBigerOneDepthCa2 = tempMin < tempDepth_one_ca2;
tempMax(markMaxBigerOneDepthCa2) = tempDepth_one_ca2( markMaxBigerOneDepthCa2);%NaN;
tempMin(markMinBigerOneDepthCa2) = tempDepth_one_ca2( markMinBigerOneDepthCa2);
tempMin(markMaxBigerOneDepthCa2) = tempDepth_zero_ca2( markMaxBigerOneDepthCa2);%NaN;
MaxD1(mark1_0) = tempMax;
MinD1(mark1_0) = tempMin;

% calculate the projection position
x1CaMax3D = inv(defaultPara.InrinsicK1)*(x1.*repmat(MaxD1,3,1)); % 3-D position in camera 1 coordinate (3 by n)
x1CaMin3D = inv(defaultPara.InrinsicK1)*(x1.*repmat(MinD1,3,1)); % 3-D position in camera 1 coordinate (3 by n)
x1CaMaxHomo = [ x1CaMax3D; ones(1,K1)]; % into homogenous coordinate (4 by n)
x1CaMinHomo = [ x1CaMin3D; ones(1,K1)]; % into homogenous coordinate (4 by n)
x1_2Max3D = [R1_2 T1_2]*x1CaMaxHomo; % 3-D position in camera 2 coordinate (3 by n)
x1_2MaxHomo = defaultPara.InrinsicK2*x1_2Max3D; % image homo coordinate in camera2 (3 by n)
x1_2Max = [ x1_2MaxHomo(1,:)./x1_2MaxHomo(3,:); x1_2MaxHomo(2,:)./x1_2MaxHomo(3,:)]; % image coordinate (2 by n)
x1_2Min3D = [R1_2 T1_2]*x1CaMinHomo; % 3-D position in camera 2 coordinate (3 by n)
x1_2MinHomo = defaultPara.InrinsicK2*x1_2Min3D; % image homo coordinate in camera2 (3 by n)
x1_2Min = [ x1_2MinHomo(1,:)./x1_2MinHomo(3,:); x1_2MinHomo(2,:)./x1_2MinHomo(3,:)]; % image coordinate (2 by n)
% x1_2Max(mark0_1) = x1_2Max(mark0_1) + (x1_2Max(mark0_1) - x1_2Min(mark0_1))*NegativeDepthTolerence;%Min529
% x1_2Min(mark1_0) = x1_2Min(mark1_0) + (x1_2Min(mark1_0) - x1_2Max(mark1_0))*NegativeDepthTolerence;%Min529
x1_2Max = x1_2Max + (x1_2Max - x1_2Min)*NegativeDepthTolerence;%Min529
x1_2Min = x1_2Min + (x1_2Min - x1_2Max)*NegativeDepthTolerence;%Min529

% Define Constrain (simple rectangle)
if FlagRotate
	tRAN = x1_2Max - x1_2Min;
	Ptr = tRAN(1,:) < 0;
	%theta_z = -atan(tRAN(2,:)./tRAN(1,:));
	theta_z = atan(tRAN(2,:)./tRAN(1,:));
	theta_z(Ptr) = theta_z(Ptr)+pi;
	Rc1_2 = [cos(-theta_z); -sin(-theta_z); sin(-theta_z); cos(-theta_z)];
	ConS1_2(1:2,:) = x1_2Min;
	ConS1_2(3,:) = sum( Rc1_2(1:2,:).*tRAN,1);% also may check sum( Rc1_2(3:4,:).*tRAN,1) close up to zeros to verify correct or not
%     Re1_2 = sum( Rc1_2(3:4,:).*tRAN,1)
	
	ConS1_2(4,:) = defaultPara.VertVar*max(H,V);
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
%	tempConS1_2([1 3 5 7],:) = min(max(tempConS1_2([1 3 5 7],:), 1), H); 
%	tempConS1_2([1 3 5 7]+1,:) = min(max(tempConS1_2([1 3 5 7]+1,:), 1), V); 
else
	Rc1_2 =[];
	ConS1_2(1,:) = min([ x1_2Max(1,:)+defaultPara.VertVar*H; x1_2Max(1,:)-defaultPara.VertVar*H;...
        	             x1_2Min(1,:)+defaultPara.VertVar*H; x1_2Min(1,:)-defaultPara.VertVar*H ],[],1);
	ConS1_2(1,:) = min( max( round(ConS1_2(1,:)), 1), H);
	ConS1_2(2,:) = max([ x1_2Max(1,:)+defaultPara.VertVar*H; x1_2Max(1,:)-defaultPara.VertVar*H;...
        	             x1_2Min(1,:)+defaultPara.VertVar*H; x1_2Min(1,:)-defaultPara.VertVar*H ],[],1);
	ConS1_2(2,:) = min( max( round(ConS1_2(2,:)), 1), H);
	ConS1_2(3,:) = min([ x1_2Max(2,:)+defaultPara.VertVar*V; x1_2Max(2,:)-defaultPara.VertVar*V;...
        	             x1_2Min(2,:)+defaultPara.VertVar*V; x1_2Min(2,:)-defaultPara.VertVar*V ],[],1);
	ConS1_2(3,:) = min( max( round(ConS1_2(3,:)), 1), V);
	ConS1_2(4,:) = max([ x1_2Max(2,:)+defaultPara.VertVar*V; x1_2Max(2,:)-defaultPara.VertVar*V;...
        	             x1_2Min(2,:)+defaultPara.VertVar*V; x1_2Min(2,:)-defaultPara.VertVar*V ],[],1);
	ConS1_2(4,:) = min( max( round(ConS1_2(4,:)), 1), V);
end

% ===========================================================

[V H] = size(I1);
% I2 project on I1
MaxD2 = D2*MaxRatio;
MinD2 = D2*MinRatio;
% calculate the prespective depth gives 0 prejective depth in camera 2
x2_1R_z = R2_1(3,:)*inv(defaultPara.InrinsicK2)*(x2);% it's possible that x1_2R_z is zero
mark =  x2_1R_z ~=0;
Depth_zero_ca2 = zeros(1, length(x2_1R_z));
Depth_one_ca2 = zeros(1, length(x2_1R_z));
Depth_zero_ca2(mark) = T2_1(3)./x2_1R_z(mark);
Depth_one_ca2(mark) = (1-T2_1(3))./x2_1R_z(mark);

% use zeros depth constraon to restrict the searching area more
mark0_1 = Depth_zero_ca2 >= Depth_one_ca2;
mark1_0 = ~mark0_1;

tempDepth_one_ca2 = Depth_one_ca2(mark0_1);
tempDepth_zero_ca2 = Depth_zero_ca2(mark0_1);
tempMax = MaxD2(mark0_1);
tempMin = MinD2(mark0_1);
markMaxBigerOneDepthCa2 = tempMax > tempDepth_one_ca2;
markMinBigerOneDepthCa2 = tempMin > tempDepth_one_ca2; % markMinBigerOneDepthCa2 : is one kind of outlier
tempMax(markMaxBigerOneDepthCa2) = tempDepth_one_ca2(markMaxBigerOneDepthCa2);
tempMax(markMinBigerOneDepthCa2) = tempDepth_zero_ca2(markMinBigerOneDepthCa2);%NaN; % indicator of outliers
tempMin(markMinBigerOneDepthCa2) = tempDepth_one_ca2(markMinBigerOneDepthCa2);%NaN; % indicator of outliers
MaxD2(mark0_1) = tempMax;
MinD2(mark0_1) = tempMin;
MaxD2(mark) = min(D2(mark)*MaxRatio,defaultPara.FarestDist);%Min529 new
MinD2(mark) = max(D2(mark)*MinRatio,defaultPara.Closestdist);%Min529 new

tempDepth_one_ca2 = Depth_one_ca2(mark1_0);
tempDepth_zero_ca2 = Depth_zero_ca2(mark1_0);
tempMax = MaxD2(mark1_0);
tempMin = MinD2(mark1_0);
markMaxBigerOneDepthCa2 = tempMax < tempDepth_one_ca2; % markMaxBigerOneDepthCa2 : is one kind of outlier
markMinBigerOneDepthCa2 = tempMin < tempDepth_one_ca2;
tempMax(markMaxBigerOneDepthCa2) = tempDepth_one_ca2( markMaxBigerOneDepthCa2);%NaN;
tempMin(markMinBigerOneDepthCa2) = tempDepth_one_ca2( markMinBigerOneDepthCa2);
tempMin(markMaxBigerOneDepthCa2) = tempDepth_zero_ca2( markMaxBigerOneDepthCa2);%NaN;
MaxD2(mark1_0) = tempMax;
MinD2(mark1_0) = tempMin;
if any( MaxD2< MinD2)
   disp('error') 
end    
% calculate the projection position
x2CaMax3D = inv(defaultPara.InrinsicK2)*(x2.*repmat(MaxD2,3,1)); % 3-D position in camera 2 coordinate (3 by n)
x2CaMin3D = inv(defaultPara.InrinsicK2)*(x2.*repmat(MinD2,3,1)); % 3-D position in camera 2 coordinate (3 by n)
x2CaMaxHomo = [ x2CaMax3D; ones(1,K2)]; % into homogenous coordinate (4 by n)
x2CaMinHomo = [ x2CaMin3D; ones(1,K2)]; % into homogenous coordinate (4 by n)
x2_1Max3D = [R2_1 T2_1]*x2CaMaxHomo; % 3-D position in camera 1 coordinate (3 by n)
x2_1MaxHomo = defaultPara.InrinsicK1*x2_1Max3D; % image homo coordinate in camera1 (3 by n)
x2_1Max = [ x2_1MaxHomo(1,:)./x2_1MaxHomo(3,:); x2_1MaxHomo(2,:)./x2_1MaxHomo(3,:)]; % image coordinate (2 by n)
x2_1Min3D = [R2_1 T2_1]*x2CaMinHomo; % 3-D position in camera 1 coordinate (3 by n)
x2_1MinHomo = defaultPara.InrinsicK1*x2_1Min3D; % image homo coordinate in camera1 (3 by n)
x2_1Min = [ x2_1MinHomo(1,:)./x2_1MinHomo(3,:); x2_1MinHomo(2,:)./x2_1MinHomo(3,:)]; % image coordinate (2 by n)
% x2_1Max(mark0_1) = x2_1Max(mark0_1) + (x2_1Max(mark0_1) - x2_1Min(mark0_1))*NegativeDepthTolerence;%Min529
% x2_1Min(mark1_0) = x2_1Min(mark1_0) + (x2_1Min(mark1_0) - x2_1Max(mark1_0))*NegativeDepthTolerence;%Min529
x2_1Max = x2_1Max + (x2_1Max - x2_1Min)*NegativeDepthTolerence;%Min529
x2_1Min = x2_1Min + (x2_1Min - x2_1Max)*NegativeDepthTolerence;%Min529
% Define Constrain (simple rectangle)
if FlagRotate
	tRAN = x2_1Max - x2_1Min;
	Ptr = tRAN(1,:) < 0;
	%theta_z = -atan(tRAN(2,:)./tRAN(1,:));
	theta_z = atan(tRAN(2,:)./tRAN(1,:));
	theta_z(Ptr) = theta_z(Ptr)+pi;
	Rc2_1 = [cos(-theta_z); -sin(-theta_z); sin(-theta_z); cos(-theta_z)];
	ConS2_1(1:2,:) = x2_1Min;
	ConS2_1(3,:) = [sum( Rc2_1(1:2,:).*tRAN,1)];
%     Re2_1 = sum( Rc2_1(3:4,:).*tRAN,1)
	
	ConS2_1(4,:) = defaultPara.VertVar*max(H,V);
        % [x_origon; y_origin; x_bound, y_bound(defaultPara.VertVar*max(H,V))]
	
        % generate plot point
	tempConS2_1 = [ ConS2_1(3,:); ConS2_1(4,:); ...
                        zeros(1,size(ConS2_1,2)); ConS2_1(4,:);...
			zeros(1,size(ConS2_1,2)); -ConS2_1(4,:);...
			ConS2_1(3,:); -ConS2_1(4,:)];
% 	Rc2_1_transpose = [cos(theta_z); -sin(theta_z); sin(theta_z); cos(theta_z)];
    Rc2_1_transpose = [cos(theta_z); -sin(theta_z); sin(theta_z); cos(theta_z)];
	tempConS2_1(1:2,:) = [sum( Rc2_1_transpose(1:2,:).*tempConS2_1(1:2,:), 1); ...
			      sum( Rc2_1_transpose(3:4,:).*tempConS2_1(1:2,:), 1)];
	tempConS2_1(3:4,:) = [sum( Rc2_1_transpose(1:2,:).*tempConS2_1(3:4,:), 1); ...
			      sum( Rc2_1_transpose(3:4,:).*tempConS2_1(3:4,:), 1)];
	tempConS2_1(5:6,:) = [sum( Rc2_1_transpose(1:2,:).*tempConS2_1(5:6,:), 1); ...
			      sum( Rc2_1_transpose(3:4,:).*tempConS2_1(5:6,:), 1)];
	tempConS2_1(7:8,:) = [sum( Rc2_1_transpose(1:2,:).*tempConS2_1(7:8,:), 1); ...
			      sum( Rc2_1_transpose(3:4,:).*tempConS2_1(7:8,:), 1)];
	tempConS2_1 = tempConS2_1 + repmat( x2_1Min, 4, 1);
    RoughConS2_1 = [ min( tempConS2_1([ 1 3 5 7],:), [], 1); max( tempConS2_1([ 1 3 5 7],:), [], 1);...
                    min( tempConS2_1([ 2 4 6 8],:), [], 1); max( tempConS2_1([ 2 4 6 8],:), [], 1)];
%	tempConS2_1([1 3 5 7],:) = min(max(tempConS2_1([1 3 5 7],:), 1), H); 
%	tempConS2_1([1 3 5 7]+1,:) = min(max(tempConS2_1([1 3 5 7]+1,:), 1), V); 
else
	Rc2_1 = [];
	ConS2_1(1,:) = min([ x2_1Max(1,:)+defaultPara.VertVar*H; x2_1Max(1,:)-defaultPara.VertVar*H;...
        	             x2_1Min(1,:)+defaultPara.VertVar*H; x2_1Min(1,:)-defaultPara.VertVar*H ],[],1);
	ConS2_1(1,:) = min( max( round(ConS2_1(1,:)), 1), H);
	ConS2_1(2,:) = max([ x2_1Max(1,:)+defaultPara.VertVar*H; x2_1Max(1,:)-defaultPara.VertVar*H;...
        	             x2_1Min(1,:)+defaultPara.VertVar*H; x2_1Min(1,:)-defaultPara.VertVar*H ],[],1);
	ConS2_1(2,:) = min( max( round(ConS2_1(2,:)), 1), H);
	ConS2_1(3,:) = min([ x2_1Max(2,:)+defaultPara.VertVar*V; x2_1Max(2,:)-defaultPara.VertVar*V;...
        	             x2_1Min(2,:)+defaultPara.VertVar*V; x2_1Min(2,:)-defaultPara.VertVar*V ],[],1);
	ConS2_1(3,:) = min( max( round(ConS2_1(3,:)), 1), V);
	ConS2_1(4,:) = max([ x2_1Max(2,:)+defaultPara.VertVar*V; x2_1Max(2,:)-defaultPara.VertVar*V;...
        	             x2_1Min(2,:)+defaultPara.VertVar*V; x2_1Min(2,:)-defaultPara.VertVar*V ],[],1);
	ConS2_1(4,:) = min( max( round(ConS2_1(4,:)), 1), V);
end

% ==========================================================
if FlagDisp
	T1_2_hat = [[0 -T1_2(3) T1_2(2)];...
                    [T1_2(3) 0 -T1_2(1)];...
               	    [-T1_2(2) T1_2(1) 0]];
	F = inv(defaultPara.InrinsicK2)'*T1_2_hat*R1_2*inv(defaultPara.InrinsicK1);
	if FlagRotate
		figure;
		dispMatchSearchRegin(I1, I2, x1, x2, tempConS1_2, tempConS2_1, F, ...
            x1_2Max, MaxD1, x1_2Min, MinD1, ...
            x2_1Max, MaxD2, x2_1Min, MinD2, ...
            FlagRotate, 'Stacking', 'h', 'Interactive', 0);
	else
		figure;
		dispMatchSearchRegin(I1, I2, x1, x2, ConS1_2, ConS2_1, F, FlagRotate, 'Stacking', 'v', 'Interactive', 0);
	end
end

return;
