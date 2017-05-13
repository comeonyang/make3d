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
function [pseudoDepthMap] = sortPixel(pixelReservoir,depthMapXRes,depthMapYRes)

imageXSize=2272;
imageYSize=1704;

% x_begin=min(pixelReservoir(:,1));
% x_end=max(pixelReservoir(:,1));
% y_begin=min(pixelReservoir(:,2));
% y_end=max(pixelReservoir(:,2));

% X_begin=round(min(pixelReservoir(:,4)));
% X_end=round(max(pixelReservoir(:,4)));
% Y_begin=round(min(pixelReservoir(:,5)));
% Y_end=round(max(pixelReservoir(:,5)));

% x_bin=floor((x_end-x_begin)/depthMapXRes);
% y_bin=floor((y_end-y_begin)/depthMapYRes);
x_values=round([(imageXSize/depthMapXRes):(imageXSize/depthMapXRes):imageXSize]-.5*(imageXSize/depthMapXRes));
y_values=round([(imageYSize/depthMapYRes):(imageYSize/depthMapYRes):imageYSize]-.5*(imageYSize/depthMapYRes));

for i=1:depthMapXRes
    for j=1:depthMapYRes
        [a,temp]=min(((pixelReservoir(:,4)-x_values(i)).^2)+((pixelReservoir(:,5)-y_values(j)).^2));
        x_error = (pixelReservoir(temp,4)-x_values(i));
        y_error = (pixelReservoir(temp,5)-y_values(j));
        pseudoDepthMap(i,j,:)=[pixelReservoir(temp,1:3) x_error y_error a];
    end
end    
