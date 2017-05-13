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
function [ Rc1_2, Rc2_1, ConS1_2, ConS2_1, RoughConS1_2, RoughConS2_1]=EffCalMatchSearchRegin(defaultPara, ScaleImg1, ScaleImg2, x1, x2, R, T, D1, D2, FlagDist);

% This function Calculate the Constrain of SurfMatch search space 
% with different information given:
% 1) Estimated Rotation and Translation matrix and depths
% 2) Estimated Rotation and Translation matrix (no depths)

% Return:
% Rc1_2 - (4 x length(x1) : the 4 element of each column is the vectorized rotation matrix) 
%		constrain for Img1 as Target, Img2 as Field
% Rc2_1 - (4 x length(x2) : the 4 element of each column is the vectorized rotation matrix) 
%		constrain for Img2 as Target, Img1 as Field
% ConS1_2 - (4 x length(x1)) constrain for Img1 as Target, Img2 as Field
%		ConS1_2([1 2],:) - reference corner for the constrain square (x y)
%		ConS1_2([3],:) - sqare width along the epipolar line
%		ConS1_2([4],:) - sqare height othorgonal to the epipolar line
% ConS2_1 - (4 x length(x2)) constrain for Img2 as Target, Img1 as Field
%		ConS2_1([1 2],:) - reference corner for the constrain square (x y)
%		ConS2_1([3],:) - sqare width along the epipolar line
%		ConS2_1([4],:) - sqare height othorgonal to the epipolar line
% RoughConS1_2 - (4 x length(x1)) constrain for Img1 as Target, Img2 as Field (not allow rotation)
%		RoughConS1_2([1 2 3 4],:) - [xmax; xmin; ymax; ymin];
% RoughConS2_1 - (4 x length(x2)) constrain for Img2 as Target, Img1 as Field (not allow rotation)
%		RoughConS2_1([1 2 3 4],:) - [xmax; xmin; ymax; ymin];

% initialize Parameter
NegativeDepthTolerence = defaultPara.NegativeDepthTolerence;
MaxRatio =defaultPara.MaxRatio;%300;
MinRatio = 1/MaxRatio;
HeightImg1 = defaultPara.VertVar*max( ScaleImg1);
HeightImg2 = defaultPara.VertVar*max( ScaleImg2);

K1 = size(x1,2);
K2 = size(x2,2);
x1 = [x1; ones(1,K1)];
x2 = [x2; ones(1,K2)];
x1_calib = inv(defaultPara.InrinsicK1)*x1;
x2_calib = inv(defaultPara.InrinsicK2)*x2;

if isempty( R) || isempty( T)
	Rc1_2 = zeros(4, K1);
	Rc1_2( [1 4],:) = 1;
	Rc2_1 = zeros(4, K2);
	Rc2_1( [1 4],:) = 1;
	ConS1_2(1,:) = zeros(1, K1);
	ConS1_2(2,:) = ScaleImg2(2)/2*ones(1, K1);
	ConS1_2(3,:) = ScaleImg2(1);
	ConS1_2(4,:) = ScaleImg2(2)/2;
	ConS2_1(1,:) = zeros(1, K2);
	ConS2_1(2,:) = ScaleImg1(2)/2*ones(1, K2);
	ConS2_1(3,:) = ScaleImg1(1);
	ConS2_1(4,:) = ScaleImg1(2)/2;
	RoughConS1_2(1,:) = ScaleImg2(1)*ones(1,K1);
	RoughConS1_2(2,:) = 0;
	RoughConS1_2(3,:) = ScaleImg2(2)*ones(1,K1);
	RoughConS1_2(4,:) = 0;
	RoughConS2_1(1,:) = ScaleImg1(1)*ones(1,K2);
	RoughConS2_1(2,:) = 0;
	RoughConS2_1(3,:) = ScaleImg1(2)*ones(1,K2);
	RoughConS2_1(4,:) = 0;
else
	R1_2 = R(1:3,:);
	R2_1 = R(4:6,:);
	T1_2 = T(1:3,:);
	T2_1 = T(4:6,:);

	T1_2_hat = [[0 -T1_2(3) T1_2(2)];...
			    [T1_2(3) 0 -T1_2(1)];...
			    [-T1_2(2) T1_2(1) 0]];
	E = T1_2_hat*R1_2;
	F = inv(defaultPara.InrinsicK2')*E*inv(defaultPara.InrinsicK1)

	% I1 project on I2 ==========================================
	% 1) calculated the closed Depth and the Farest Depth that can be seen from Img2
	%    find Two End points of Epipolar line on Img2
	[ EndPointsImg2 ] = EndPointsFromF(F, x1, ScaleImg2);
	[ EndPointsDepthImg1(1,:) dump Error] = triangulation( defaultPara, R1_2, T1_2, [x1; inv(defaultPara.InrinsicK1)*[EndPointsImg2(1:2,:); ones(1,K1)]]);
	[ EndPointsDepthImg1(2,:) dump Error] = triangulation( defaultPara, R1_2, T1_2, [x1; inv(defaultPara.InrinsicK2)*[EndPointsImg2(3:4,:); ones(1,K1)]]);
	EndPointsDepthImg1 = sort(EndPointsDepthImg1,1); % make the EndPointsDepthImg1 in acend order from top to bottom

	% 2) prune depth range 
	if ~isempty( D1 )
		MaxD1 = D1*MaxRatio;
		MinD1 = D1*MinRatio;
		%    prune by EndPointsDepth
		MaxD1 = min(MaxD1, EndPointsDepthImg1(2,:));
		MaxD1 = max(MaxD1, EndPointsDepthImg1(1,:));
		MinD1 = max(MinD1, EndPointsDepthImg1(1,:));
		MinD1 = min(MinD1, EndPointsDepthImg1(2,:));
	else
		MaxD1 = EndPointsDepthImg1(2,:);
		MinD1 = EndPointsDepthImg1(1,:);
	end
	%   prune by additional constrain ========OPtional
	MaxD1 = min(MaxD1, defaultPara.FarestDist);
	MinD1 = max(MinD1, defaultPara.Closestdist); 
	%   ==============================================

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

	% expand the search space a little bit in case the R and T are not accurate enough
	x1_2Max = x1_2Max + (x1_2Max - x1_2Min)*NegativeDepthTolerence;%Min529
	x1_2Min = x1_2Min + (x1_2Min - x1_2Max)*NegativeDepthTolerence;%Min529

	% Define Constrain (simple rectangle)
	[ Rc1_2, ConS1_2, RoughConS1_2 ]=Points2SqareConstrain( [ x1_2Max; x1_2Min], HeightImg1);

	% ===========================================================

	% I2 project on I1 ==========================================
	% 1) calculated the closed Depth and the Farest Depth that can be seen from Img1
	%    find Two End points of Epipolar line on Img2
	[ EndPointsImg1 ] = EndPointsFromF(F', x2, ScaleImg1);
	[ EndPointsDepthImg2(1,:) dump Error] = triangulation( defaultPara, R2_1, T2_1, [x2; inv(defaultPara.InrinsicK1)*[EndPointsImg1(1:2,:); ones(1,K2)]]);
	[ EndPointsDepthImg2(2,:) dump Error] = triangulation( defaultPara, R2_1, T2_1, [x2; inv(defaultPara.InrinsicK2)*[EndPointsImg1(3:4,:); ones(1,K2)]]);
	EndPointsDepthImg2 = sort(EndPointsDepthImg2,1); % make the EndPointsDepthImg1 in acend order from top to bottom

	% 2) prune depth range 
	if ~isempty( D2)
		MaxD2 = D2*MaxRatio;
		MinD2 = D2*MinRatio;
		%    prune by EndPointsDepth
		MaxD2 = min(MaxD2, EndPointsDepthImg2(2,:));
		MaxD2 = max(MaxD2, EndPointsDepthImg2(1,:));
		MinD2 = max(MinD2, EndPointsDepthImg2(1,:));
		MinD2 = min(MinD2, EndPointsDepthImg2(2,:));
	else
		MaxD2 = EndPointsDepthImg2(2,:);
		MinD2 = EndPointsDepthImg2(1,:);
	end
	%   prune by additional constrain ========OPtional
	MaxD2 = min(MaxD2, defaultPara.FarestDist);
	MinD2 = max(MinD2, defaultPara.Closestdist); 
	%   ==============================================

	% calculate the projection position
	x2CaMax3D = inv(defaultPara.InrinsicK2)*(x2.*repmat(MaxD2,3,1)); % 3-D position in camera 1 coordinate (3 by n)
	x2CaMin3D = inv(defaultPara.InrinsicK2)*(x2.*repmat(MinD2,3,1)); % 3-D position in camera 1 coordinate (3 by n)
	x2CaMaxHomo = [ x2CaMax3D; ones(1,K2)]; % into homogenous coordinate (4 by n)
	x2CaMinHomo = [ x2CaMin3D; ones(1,K2)]; % into homogenous coordinate (4 by n)
	x2_1Max3D = [R2_1 T2_1]*x2CaMaxHomo; % 3-D position in camera 2 coordinate (3 by n)
	x2_1MaxHomo = defaultPara.InrinsicK2*x2_1Max3D; % image homo coordinate in camera2 (3 by n)
	x2_1Max = [ x2_1MaxHomo(1,:)./x2_1MaxHomo(3,:); x2_1MaxHomo(2,:)./x2_1MaxHomo(3,:)]; % image coordinate (2 by n)
	x2_1Min3D = [R2_1 T2_1]*x2CaMinHomo; % 3-D position in camera 2 coordinate (3 by n)
	x2_1MinHomo = defaultPara.InrinsicK2*x2_1Min3D; % image homo coordinate in camera2 (3 by n)
	x2_1Min = [ x2_1MinHomo(1,:)./x2_1MinHomo(3,:); x2_1MinHomo(2,:)./x2_1MinHomo(3,:)]; % image coordinate (2 by n)

	% expand the search space a little bit in case the R and T are not accurate enough
	x2_1Max = x2_1Max + (x2_1Max - x2_1Min)*NegativeDepthTolerence;%Min529
	x2_1Min = x2_1Min + (x2_1Min - x2_1Max)*NegativeDepthTolerence;%Min529

	% Define Constrain (simple rectangle)
	[ Rc2_1, ConS2_1, RoughConS2_1 ]=Points2SqareConstrain( [ x2_1Max; x2_1Min], HeightImg2);

	% ===========================================================
	if FlagDisp
			%figure;
			%dispMatchSearchRegin(I1, I2, x1, x2, tempConS1_2, tempConS2_1, F, ...
		    %x1_2Max, MaxD1, x1_2Min, MinD1, ...
		    %x2_1Max, MaxD2, x2_1Min, MinD2, ...
		    %FlagRotate, 'Stacking', 'h', 'Interactive', 0);
			figure;
			dispMatchSearchRegin(I1, I2, x1, x2, tempConSConS1_2, tempConSConS2_1, F, FlagRotate, 'Stacking', 'v', 'Interactive', 0);
		end
	end
end

% return
