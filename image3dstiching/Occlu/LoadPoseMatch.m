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
function [Pair Matches, fail]=LoadPoseMatch(Fdir, ImgName1, ImgName2)

% This function load the xxx_xxx_PoseMatch.mat proberly by tring the two order automatically
[status1, result1] = system(['ls ' Fdir '/data/' ImgName1 '_' ImgName2 '_PoseMatch.mat']);
[status2, result2] =  system(['ls ' Fdir '/data/' ImgName2 '_' ImgName1 '_PoseMatch.mat']);
if ~status1
    load([Fdir '/data/' ImgName1 '_' ImgName2 '_PoseMatch.mat']);
	% Cool then just load the .mat file
elseif ~status2
    load([Fdir '/data/' ImgName2 '_' ImgName1 '_PoseMatch.mat']);
	% Well we just need ot switch the variable    
	Pair_tmp = Pair;
	Pair.R = Pair_tmp.R';
	Pair.T = -Pair_tmp.R'*Pair_tmp.T;
	Pair.Xim = Pair_tmp.Xim([3 4 1 2],:);
	Pair.DepthScale = Pair_tmp.DepthScale([2 1],:);
	Pair.lamda = Pair_tmp.lamda([2 1],:);
	Pair.t = Pair_tmp.t;
	
	Matches = Matches([2 1],:);
	
else % no direct match between ImgName1 and ImgName2
     % so change to find any match including ImgName1 and ImgName2 seperately
	% for ImgName1
	% No used anymore =======================================================================
     if false
	file1 = dir([Fdir '/data/']);
	for i = 1:length(file1)
		if ~isempty( strfind(file1(i).name,ImgName1))

			load([Fdir '/data/' file1(i).name]);
			if strfind(file1(i).name,ImgName1) == 1
				Pair_tmp.Xim1 =	Pair.Xim(1:2,:);
				Pair_tmp.DepthScale1 = Pair.DepthScale(1);
				Pair_tmp.lamda1 = Pair.lamda(1,:);
			else
				Pair_tmp.Xim1 =	Pair.Xim(3:4,:);
				Pair_tmp.DepthScale1 = Pair.DepthScale(1);
				Pair_tmp.lamda1 = Pair.lamda(1,:);
			end
		end
	end
     end
	% ===================================================================================

     % just simply use empty matches;
	Pair.Xim = zeros(4,0);
	Pair.DepthScale = zeros(0,2);
	Pair.lamda = zeros(2,0);	

	Matches = zeros(2,0);
	fail = 0;		
	disp('No Match been tested');
end

return;
