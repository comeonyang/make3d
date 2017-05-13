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
function [depthBins,countBuf,count1,count2,count3,count4]=trimDepthBuf(depthBins,countBuf,panBins,sweepBins)

%% first trim columns

count4=0;
for j=1:sweepBins
    j=1;
    tmp=sum(depthBins(:,j)');
    if(tmp==0)
        depthBins=depthBins(:,[j+1:end]);
        countBuf=countBuf(:,[j+1:end]);
        count4=count4+1;
    else
        break;
    end
end
count3=0;
for j=size(depthBins,2):-1:1
    j=size(depthBins,2);
    tmp=sum(depthBins(:,j)');
    if(tmp==0)
        depthBins=depthBins(:,1:end-1);
        countBuf=countBuf(:,1:end-1);
        count3=count3+1;
    else
        break;
    end
end

%% then trim rows
count2=0;
for i=1:panBins
    i=1;
    tmp=sum(depthBins(i,:)');
    if(tmp==0)
        depthBins=depthBins([i+1:end],:);
        countBuf=countBuf([i+1:end],:);
        count2=count2+1;
    else
        break;
    end
end
count1=0;
for i=size(depthBins,1):-1:1
    i=size(depthBins,1);
    tmp=sum(depthBins(i,:)');
    if(tmp==0)
        depthBins=depthBins([1:end-1],:);
        countBuf=countBuf([1:end-1],:);
        count1=count1+1;
    else
        break;
    end
end
