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
function [depthBins,countBuf]=checkIfValid(depthBins,countBuf,panBins,sweepBins,minDepth)

for i=1:panBins
    for j=1:sweepBins
        if ((countBuf(i,j)==0)||(depthBins(i,j))==0) %% take the nearest non-zero value in the same column
            check_col=j;
            [temp1,temp2]=find(depthBins(:,check_col)');
%             while (length(temp1)==0)
%                 if (check_col<sweepBins)
%                     check_col=check_col+1;
%                 else
%                     check_col=check_col-1;
%                 end
%                 [temp1,temp2]=find(depthBins(:,check_col)');                
%             end            
            [temp3,temp4]=min(abs(temp2-i));
% 			if ((min(depthBins(temp2,check_col)./countBuf(temp2,check_col))<minDepth)||(countBuf(temp2(temp4),check_col)<1))
%                 waithere=1;
% 			end
            depthBins(i,j)=depthBins(temp2(temp4),check_col);
            countBuf(i,j)=countBuf(temp2(temp4),check_col);
        end 
    end
end
depthBins=depthBins./countBuf;
