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
function [] = generate_vrml_faceset_triangle(vrmlDirectory, vrmlFilename, imgDirectory, imageFilename, leftTriangle, rightTriangle)
% this function play with the VRML
% using only FaceSet

% if nargin < 6
%     a = 0.70783777 %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
%     b = 0.946584169%0.085; % vertical physical size of image plane normalized to focal length (in meter)
%     Ox = -0.010727086; % camera origin offset from the image center in horizontal direction
%     Oy = -0.0111130176; % camera origin offset from the image center in vertical direction
% elseif nargin < 7
%     b = a;
%     Ox = -0.010727086; % camera origin offset from the image center in horizontal direction
%     Oy = -0.0111130176; % camera origin offset from the image center in vertical direction
% end
% 
% % define global variable
% global GeneralDataFolder ScratchDataFolder LocalFolder;
% 
% vrml_filename = ['test_faceset_triangle' filename '_' DepthDirectory '.wrl'];
% LowResImgIndexSuperpixelSepTemp = LowResImgIndexSuperpixelSep;
imgDirectoryMod=strcat('.',imgDirectory);
A = imread([imgDirectoryMod imageFilename]);

% [VertYSize HoriXSize] = size(A);
% nu_patch = VertYSize* HoriXSize;
% PlaneParameterTureTemp = PlaneParameterTure;

% calculate ray
% ray = GenerateRay(HoriXSize,VertYSize,'center'); %[ horiXSizeLowREs VertYSizeLowREs 3]

% calculate Position3D
% DepthTureProj = 1./sum(permute(reshape(PlaneParameterTureTemp(1:3,LowResImgIndexSuperpixelSepTemp),3,VertYSize,[]),[2 3 1]).*ray,3);
% Position3D = im_cr2w_cr(DepthTureProj,ray);

% calculate Position3DCoord
Position3DCoord = [];
for i = 1:min(size(leftTriangle,1),size(rightTriangle,1))
    Position3DCoord = [Position3DCoord assignTriangleValues(leftTriangle(i,3,3:5))' assignTriangleValues(leftTriangle(i,1,3:5))' assignTriangleValues(leftTriangle(i,2,3:5))' ...
        assignTriangleValues(rightTriangle(i,2,3:5))' assignTriangleValues(rightTriangle(i,1,3:5))' assignTriangleValues(rightTriangle(i,3,3:5))'];
end

%     for j = 1:HoriXSize-1
%         Position3DCoord = [Position3DCoord Position3D(:,i,j) Position3D(:,i+1,j) Position3D(:,i+1,j+1) Position3D(:,i,j) Position3D(:,i+1,j+1) Position3D(:,i,j+1)];
%     end
% end    
Position3DCoord(3,:) = -Position3DCoord(3,:); % important to make z direction negative

% calculate PositionTexCoord
% PositionTex = permute(ray(:,:,1:2)./repmat(cat(3,a,b),[VertYSize HoriXSize 1])+0.5,[3 1 2]);
PositionTexCoord = [];
for i = 1:min(size(leftTriangle,1),size(rightTriangle,1))
    PositionTexCoord = [PositionTexCoord assignTriangleValues(leftTriangle(i,3,1:2))' assignTriangleValues(leftTriangle(i,1,1:2))' assignTriangleValues(leftTriangle(i,2,1:2))' ...
        assignTriangleValues(rightTriangle(i,2,1:2))' assignTriangleValues(rightTriangle(i,1,1:2))' assignTriangleValues(rightTriangle(i,3,1:2))'];
end
% PositionTexCoord = PositionTexCoord./([2272 1704]');
PositionTexCoord = PositionTexCoord./(repmat([2272 1704]',1,length(PositionTexCoord)));
pathName = pwd;
imageForVRMLfile=strcat(pathName,imgDirectory,imageFilename);

    % inital header
disp('writing vrml..');
fp = fopen([vrmlDirectory vrmlFilename],'w');

fprintf(fp, '#VRML V2.0 utf8\n');

% add navigate_info
fprintf(fp, 'NavigationInfo {\n');
fprintf(fp, '  headlight TRUE\n');
fprintf(fp, '  type ["FLY", "ANY"]}\n\n');

% add viewpoint
fprintf(fp, 'Viewpoint {\n');
fprintf(fp, '    position        0 0.0 0.0\n');
fprintf(fp, '    orientation     0 0 0 0\n');
fprintf(fp, '    fieldOfView     0.1\n');
fprintf(fp, '    description "Original"}\n');

% add Shape for texture faceset
fprintf(fp, 'Shape{\n');
fprintf(fp, '  appearance Appearance {\n');
fprintf(fp, ['   texture ImageTexture { url "' imageForVRMLfile '" }\n']);
fprintf(fp, '  }\n');
fprintf(fp, '  geometry IndexedFaceSet {\n');
fprintf(fp, '    coord Coordinate {\n');

% insert coordinate in 3d
% =======================
fprintf(fp, '      point [ \n');
fprintf(fp, '        %f %f %f,\n',Position3DCoord);
fprintf(fp, '      ]\n');
fprintf(fp, '    }\n');

% insert coordinate index in 3d
fprintf(fp, '    coordIndex [\n');
Position3DCoordOrdeerIndex = 0:(size(Position3DCoord,2)-1);
fprintf(fp, '              %g %g %g -1,\n',Position3DCoordOrdeerIndex);
fprintf(fp, '    ]\n');

% insert texture coordinate
fprintf(fp, '	 texCoord TextureCoordinate {\n');
fprintf(fp, '	   point [\n');
fprintf(fp, '              %g %g,\n',PositionTexCoord);
fprintf(fp, '	     ]\n');
fprintf(fp, '	 }\n');
% ==================================
fprintf(fp, '  }\n');
fprintf(fp, '}\n');

% close the file
fclose(fp);
disp('done :)');
