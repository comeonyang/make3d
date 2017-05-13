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
function [f1, f2, matches] = readSurfMatches(file1, file2, dir, Type, DenseFlag, Flag128, NumCol)
% This function reads in surf features and matches for a pair of images.
% file1/2 contain the features,  contians the matches.
% if dir is passed, find the featues at dir/surf/file1.surf
% and matches at dir/surf_matches/file1-file2.match
if nargin < 5
    DenseFlag = false;
    Flag128 = true;
    NumCol = 2;
elseif nargin < 6
    Flag128 = true;
    NumCol = 2;
elseif nargin <7
    NumCol = 2;
end

FlagSwitch = 0;
% initialize surf1 surf2 and match
if DenseFlag
    surf1 = [dir '/surf/' file1 '.surfDense'];
    surf2 = [dir '/surf/' file2 '.surfDense'];
    match = [dir '/surf_matches/' file1 '-' file2 '.match' Type];
    f = fopen(match);
    if f==-1
        match = [dir '/surf_matches/' file2 '-' file1 '.match' Type];
        FlagSwitch = 1;
        f = fopen(match);
    end
else
    if Flag128
        surf1 = [dir '/surf/' file1 '.surf'];
        surf2 = [dir '/surf/' file2 '.surf'];
        match = [dir '/surf_matches/' file1 '-' file2 '.match' Type];
    else
        surf1 = [dir '/surf/' file1 '.surf_64'];
        surf2 = [dir '/surf/' file2 '.surf_64'];
        match = [dir '/surf_matches/' file1 '-' file2 '.match_64'];
    end
    f = fopen(match);
    if f==-1
        FlagSwitch = 1;
        if Flag128
            surf2 = [dir '/surf/' file2 '.surf'];
            surf1 = [dir '/surf/' file1 '.surf'];
            match = [dir '/surf_matches/' file2 '-' file1 '.match' Type];
        else
            surf2 = [dir '/surf/' file2 '.surf_64'];
            surf1 = [dir '/surf/' file1 '.surf_64'];
            match = [dir '/surf_matches/' file2 '-' file1 '.match_64'];
        end
        f = fopen(match);
    end
end

if f == -1
	f1 = [];
	f2 = [];
	matches = [];
    return;
end

matches = fscanf(f, '%g', [NumCol, inf]); %2xNumMatches
matches = matches+1;
fclose(f);
if FlagSwitch
   matches = flipud(matches);
end


% read in first set of features
f = fopen(surf1);
g = fgetl(f); %skip first 2 lines
features = str2num(g)+5; % dynamic deciding features length from data
g = fgetl(f);
%70xnumPoints, xy positions are points(1:2, :) for SURF-64
%134xnumPoints, xy positions are points(1:2, :) for SURF-128
frames1 = fscanf(f, '%f', [features, inf]); 
fclose(f);

% read in second set of features
f = fopen(surf2);
g = fgetl(f); %skip first 2 lines
g = fgetl(f);
frames2 = fscanf(f, '%f', [features, inf]); %70xnumPoints, xy positions are points(1:2, :)
fclose(f);

f1 = frames1(1:2, :);
f2 = frames2(1:2, :);
