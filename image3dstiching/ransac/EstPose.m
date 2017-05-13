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
function [ R, T, lamda1, lamda2, inlier, Error] = EstPose( defaultPara, E, x, PriorDepthRatio, PriorR)

% This function estimate the Pose ( R rotation, T translation) from E
% After finding R and T, it also triangulate the depth(lamda) when 
% T is unit length
% Input:
%         E - essential matrix (for calibrated camera)
%         x - calibrated matches point in both images
%
% Return
%         R - rotation matrix
%         T - translation unit vector
%         lamda1/2 - triangulated depth for image 1 and 2 when T is unit length

LowThre  = 0.2;
 [U S V] =svd(E);
%  pause
 Rz_pos = [ [0 -1 0];...
            [1 0 0];...
            [0 0 1]];

 R1 = U*Rz_pos'*V';
 T1 = U(:,end);
 
 R2 = U*Rz_pos*V';
 T2 = -T1;

 % eight/four possible combination to be the solution.
 % check positive depth constrain
 lamda1 = [];
 lamda2 = [];
 R = [];
 T = [];
 % Building RPool
 count = 1;
 RPoolAll = {R1,R2,-R1,-R2};
 for j = 1:length(RPoolAll)
	% test UpSide Dwon
	NewZ = RPoolAll{j}*[0; 1; 0];
	if NewZ(2) > 0
		if ~DetectImproperRotation(RPoolAll{j}) % test not improperRotation
			RPool{count} = RPoolAll{j};
			count = count + 1;
		end
	end
 end
 TPool = {T1,T2};

 count = 1;
 for j = 1:length(RPool)
	for k = 1:length(TPool)
		 [lamda1Candidate{count}, lamda2Candidate{count} Error{count}] = triangulation( defaultPara, RPool{j}, TPool{k}, x);
		 PositiveRatio(count) = sum(lamda1Candidate{count} >0 & lamda2Candidate{count} > 0) / size(lamda1Candidate{count}, 2);
		 inlierM{count} = find(lamda1Candidate{count} >	0 & lamda2Candidate{count} > 0);
         count = count +1;
	end
 end

 % filter out some obvious bad matches
 IND = find( PositiveRatio > LowThre); 
 if length(IND) ~= 1 % multiple choice cases
	 if isempty(PriorR)

		if ~isempty(PriorDepthRatio)
			% choose the lamdaRatio close to PriorDepthRatio
			
        else
			[C I] = max( PositiveRatio( IND));
            TempIND = find( PositiveRatio( IND) == C);
            if length(TempIND) > 1
                MedainErro = median( cell2mat(Error( TempIND)'), 2);
                [Ct It] = min(MedainErro);
    			I = IND( TempIND(It));
            else
                I = IND( I);
            end    
		end	
	 else
		% choose the R closer to PriorR
		[PriorRAxis, q] = Rotation2Q(PriorR);
 		count = 1;
		for i = IND
   				[RAxis(count,:), q] = Rotation2Q(RPool{ floor( (i-1)/length(TPool))+1});
                count = count + 1;
		end
		DistRotationAxis = sqrt( sum( ( RAxis(:,1:3) - repmat(PriorRAxis(1,1:3), length(IND), 1)).^2, 2));
		[Ctemp I] = min(DistRotationAxis);
		I = IND(I);
     end 
     
 else
     I = IND; % single choice cases
 end

 if ~isempty(I)
    [i j] = ind2sub([ length(TPool) length(RPool)],I);
    R = RPool{j};
    T = TPool{i};
    inlier = inlierM{I};
    lamda1 = lamda1Candidate{I};
    lamda2 = lamda2Candidate{I};
    Error = Error{I};
 else
    R = [];
    T = [];   
    inlier = [];
    lamda1 = [];
    lamda2 = [];
    Error = [];
 end
 
 if isempty(R)
	 disp('all four solution failed');
 end
 return;
