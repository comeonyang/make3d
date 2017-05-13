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
Fdir ='/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/ImgWithLaserData/'
Img1 = 'img-manmade6-p-107t0'
Img2 = 'img-manmade7-p-343t0'

Default.fy = 2400.2091651084;
Default.fx = 2407.3312729885838;
Default.Ox = 1110.7122391785729;%2272/2; %
Default.Oy = 833.72104535435108;%1704/2; %
Default.a_default = 2272/Default.fx; %0.70783777; %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
Default.b_default = 1704/Default.fy; %0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
Default.Ox_default = 1-Default.Ox/2272;%0.489272914; % camera origin offset from the image center in horizontal direction
Default.Oy_default = 1-Default.Oy/1704;%0.488886982; % camera origin offset from the image center in vertical direction

load([ Fdir '/data/depth_sph_corr-manmade6-p-107t0.mat']);
LaserDepth1 = Position3DGrid(:,:,4);

load([ Fdir '/data/depth_sph_corr-manmade7-p-343t0.mat']);
LaserDepth2 = Position3DGrid(:,:,4);
[Default.VertYNuDepth Default.HoriXNuDepth] = size( LaserDepth1);
GoodMask = LaserDepth1 <80;

load([ Fdir '/data/'  Img1 '/Multi_TestLaserImg_img-manmade6-p-107t0_NonMono.mat' ]);
TriDepth1 = model.Constrain.Depth_modified;
[IND1] = ProjPosi2Mask( [2272 1704], [55 305], Pair.Xim([1 2],:));
[ReScale] = Rescale(TriDepth1, LaserDepth1(IND1));
TriDepth1 = TriDepth1*ReScale;
MonoDepth1 = model.Depth.FitDepth;
[ReScale] = Rescale(MonoDepth1, LaserDepth1);
MonoDepth1 = ReScale*MonoDepth1;
JointDepth1 = model.PlaneParaInfo.FitDepth;
[ReScale] = Rescale(JointDepth1, LaserDepth1);
JointDepth1 = ReScale*JointDepth1;

[ Mrms Mabs_m Mrms_f Mabs_m_f] = errorMetric( MonoDepth1(GoodMask), LaserDepth1(GoodMask))
[ Jrms Jabs_m Jrms_f Jabs_m_f] = errorMetric( JointDepth1(GoodMask), LaserDepth1(GoodMask))
[ Trms Tabs_m Trms_f Tabs_m_f] = errorMetric( TriDepth1(:), LaserDepth1(IND1)')

% [ rms abs_m rms_f abs_m_f] = errorMetric(D1, LaserDepth1)

load([ Fdir '/data/'  Img2 '/Multi_TestLaserImg_img-manmade7-p-343t0_NonMono.mat' ]);
GoodMask = LaserDepth2 <80;
TriDepth2 = model.Constrain.Depth_modified;
[IND2] = ProjPosi2Mask( [2272 1704], [55 305], Pair.Xim([3 4],:));
[ReScale] = Rescale(TriDepth2, LaserDepth2(IND2));
TriDepth2 = TriDepth2*ReScale;

MonoDepth2 = model.Depth.FitDepth;
[ReScale] = Rescale(MonoDepth2, LaserDepth2);
MonoDepth2 = ReScale*MonoDepth2;
JointDepth2 = model.PlaneParaInfo.FitDepth;
[ReScale] = Rescale( JointDepth2, LaserDepth2);
JointDepth2 = ReScale*JointDepth2;

[ Mrms Mabs_m Mrms_f Mabs_m_f] = errorMetric( MonoDepth2(GoodMask), LaserDepth2(GoodMask))
[ Jrms Jabs_m Jrms_f Jabs_m_f] = errorMetric( JointDepth2(GoodMask), LaserDepth2(GoodMask))
[ Trms Tabs_m Trms_f Tabs_m_f] = errorMetric( TriDepth2(:), LaserDepth2(IND2)')

% rescale
scaleMono1 = sdpvar(1,1);
sol = solvesdp([],norm( scaleMono1*MonoDepth1(:) - LaserDepth1(:),1));
