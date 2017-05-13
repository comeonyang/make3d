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

% ============== camera intrinsic parameter
% Px = 1704 pixels
% Py = 2272 pixels
% Sx : physical size of the CCD in x direction
% Sy : physical size of the CCD in y direction
% f: focal lenghth
Default.fy = 2400.2091651084; % (Py / Sy) * f means the length of focal in terms of pixel number in y direction
Default.fx = 2407.3312729885838; % (Px / Sx) * f means the length of focal in terms of pixel number in x direction
Default.Ox = 1110.7122391785729; % camera center of horizontal direction im pixel coordinate ( left as zero)
Default.Oy = 833.72104535435108; % camera center of vertical direction im pixel coordinate (top as zero)
%Default.Ox_default = Default.Ox/2272; % 0.4889; % camera origin offset from the image center in horizontal direction
Default.Ox_default = 1-Default.Ox/2272; % 0.4889; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;% 0.5107 % camera origin offset from the image center in vertical direction
Default.b_default = 1704/Default.fy; %0.7078 (Sx / f) horizontal physical size of CCD normalized to focal length
Default.a_default = 2272/Default.fx; %0.9466 (Sy / f) vertical physical size of CCD normalized to focal length

% K (camera intrinsic matrix): 
% [ f*Px/Sx 0 0; 0 f*Py/Sy 0; 0 0 1] for pixel coordinate
% ===================================================


Default.SegVertYSize = 900%1200;
Default.SegHoriXSize = 1200%900;
Default.VertYNuPatch = 61%55;
Default.HoriXNuPatch = 55%61;%305;%61;
Default.VertYNuDepth = 55;
Default.HoriXNuDepth = 305;
Default.PopUpHoriX = 600;  
Default.PopUpVertY = 800;
Default.batchSize = 10;
Default.NuRow_default = 55;
Default.WeiBatchSize = 5;
Default.TrainVerYSize = 1704%2272;
Default.TrainHoriXSize = 2272%1704;
Default.MempryFactor =2;


