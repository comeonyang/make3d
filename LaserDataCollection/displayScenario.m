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
function displayScenario(pixelReservoir,imgDirectory)

len=size(pixelReservoir,1);
wid=size(pixelReservoir,2);
imagePanAnglePrev=1000;
for i=1:len
    for j=1:wid
        imagePanAngle=pixelReservoir(i,j,2);
        if (imagePanAngle~=imagePanAnglePrev)
            count=1;
            %%load the image in A
            if (imagePanAnglePrev ~= 1000)
                figure;
                subplot(1,2,1)
                imagesc(A), axis equal; 
                %% also display ima
                B = rgb2ycbcr(A);
				scaling_factor = floor(256/max(max(ima)));
				B(:,:,2) = scaling_factor*ima;
                subplot(1,2,2)
                imagesc(B), axis equal; 
                clear ima;
            end
            imagePanAnglePrev=imagePanAngle;
            dirList = dir([imgDirectory strcat('*',num2str(abs(imagePanAngle)),'*.jpg')]);
            if (length(dirList)==0)
                problem_list=1
            end
			imageFilename = dirList(1).name;
            A = imread([imgDirectory imageFilename]);
            ima=zeros(size(A,1), size(A,2));
		else
            ima(floor(pixelReservoir(i,j,3)), floor(pixelReservoir(i,j,4)))= pixelReservoir(i,j,1); %% x y and z
            count=count+1;
        end
    end
end
            
            
        
