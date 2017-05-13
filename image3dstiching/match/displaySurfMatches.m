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
function displaySurfMatches(dir, file1, file2, Type, DenseFlag, Flag128)

% This function display the SurfMatches between file1 and file2
% options:
%         DenseFalg = 1 if using Dense matching
%         Flag128 = 1 if using SURF-128 descriptor

if nargin < 4
    DenseFlag = false;
    Flag128 = true;
elseif nargin < 5
    Flag128 = true;
end

% get images
I1=imreadbw([dir '/pgm/' file1 '.pgm']); % function from sift
I2=imreadbw([dir '/pgm/' file2 '.pgm']); % function from sift

[f1, f2, matches] = readSurfMatches(file1, file2,  dir, Type, DenseFlag, Flag128);
% display matches	% function in image3dstiching/match
figure;
plotmatches(I1,I2,f1, f2,matches, 'Stacking', 'v', 'Interactive', 1); %change to 1 to see all matches at once (function from sift)

