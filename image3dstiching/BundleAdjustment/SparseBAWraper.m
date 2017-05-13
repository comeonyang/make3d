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
function [NewR NewT NewX_obj NewX_im dist1 dist2]=SparseBAWraper(defaultPara, R, T, x_im, X_obj, ImgInfo, FlagDisp)

% This function is the wraper function to call eucsbademo
% which perform the sparse bundle adjustment given
% Input:
%	1) defaultPara - camera intrinsic paramter in calib.txt
%	2) Rotation adn Translation - camera viewpoint in cams.txt
%	3) X_obj - object 3D position in pts.txt
% note: for the format of calib.txt cams.txt pts.txt 
%	please see README.txt in sba-1.3/demo/
%
% Return:
%	1) 

Img1 = strrep(ImgInfo(1).ExifInfo.name,'.jpg','');
Img2 = strrep(ImgInfo(2).ExifInfo.name,'.jpg','');
NumObj = length(X_obj);

% calculate qauternion
[axis_angle qauternion] = Rotation2Q(R(1:3,:));

% Write 3 .txt files in Fdir/info/
% 1) cams.txt
fp = fopen([defaultPara.Fdir '/info/cams'],'w');
fprintf(fp, '1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000\n'); % the reference camera
fprintf(fp, '%.6g ', qauternion); % print the quaternion of the second camera
fprintf(fp, '%.6g ', T(1:3)); % print the translation of the second camera
fclose(fp);

% 2) pts.txt
fp = fopen([defaultPara.Fdir '/info/pts'],'w');
fprintf(fp, '# X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ... \n'); % the reference camera
for i = 1:NumObj
	fprintf(fp, '%.6g ', X_obj(:,i));
	fprintf(fp, '2 0 ');
	fprintf(fp, '%.6g ', x_im(1:2,i));
	fprintf(fp, '1 ');
	fprintf(fp, '%.6g ', x_im(3:4,i));
	fprintf(fp, '\n');
end
fclose(fp);

% 3) calib.txt
fp = fopen([defaultPara.Fdir '/info/calib'],'w');
fprintf(fp, '%.6g ', defaultPara.InrinsicK1(1,:));
fprintf(fp, '\n');
fprintf(fp, '%.6g ', defaultPara.InrinsicK1(2,:));
fprintf(fp, '\n');
fprintf(fp, '%.6g ', defaultPara.InrinsicK1(3,:));
fprintf(fp, '\n');
fclose(fp);

%--------------------------------------------------

% call sba to solve sfm
sbaInfoPath = [defaultPara.Fdir '/info/'];
outputName = [ Img1 '-' Img2 '.out'];
system(['../third_party/sba-1.3/demo/eucsbademo ' sbaInfoPath 'cams ' sbaInfoPath 'pts ' sbaInfoPath 'calib > ' sbaInfoPath outputName]);

% Read in the result R T and X_obj x_im
[NewR NewT NewX_im NewX_obj] = readBA([sbaInfoPath outputName]);
% system(['rm ' sbaInfoPath outputName]);

temp = defaultPara.InrinsicK1*NewX_obj;
NewX_im(1:2,:) = temp(1:2,:)./repmat( temp(3,:),2,1);
temp = defaultPara.InrinsicK2*(NewR*NewX_obj + repmat(NewT,1,length(NewX_obj)));
NewX_im(3:4,:) = temp(1:2,:)./repmat( temp(3,:),2,1);

re = NewX_im - x_im;
dist1 = sqrt(sum(re(1:2,:).^2, 1));
dist2 = sqrt(sum(re(3:4,:).^2, 1));

return;
