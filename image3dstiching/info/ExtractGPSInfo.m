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
function [ImgInfo]=ExtractGPSInfo(defaultPara,Fdir, ImgInfo)

% This function extract the GPS information
% ex: Longitude Latitude Altitude (Maybe more in future)
% Then attach the Info on the ImgInfo if applicable

% initialize file identifier
temp = dir([Fdir '/info/*.ubx']);
fid = fopen([Fdir '/info/' temp.name]);
SumAltitude = 0;

PriorTime = '0';
i = 1;
% extract GPS info from .ubx file
while i <= length(ImgInfo)
    GPSWorks = true;
	TargetTime = ImgInfo(i).ExifInfo.DateTime;
    Ptr = strfind(TargetTime,' ')+1;
    TargetTime = TargetTime(Ptr:end);
    TargetTime = strrep(TargetTime,':','');
    % finding $GPGGA  
    inline_tmp = [];
	while str2num(TargetTime) > str2num(PriorTime)
		nline = fgetl(fid);
        if nline==-1
            disp('end of file');
            break;
        end
		if strcmp( nline(1:6), '$GPGGA')
            inline_tmp = nline;
 			PriorTime = nline(8:13);
                end
    end
    nline = inline_tmp;
    % locate
    Ptr = strfind(nline,',');
    % Latitude  info in $GPGGA line
    ImgInfo(i).Geo.Latitude = str2num( (nline( (Ptr(2)+1): (Ptr(2)+2) ) ) ) + str2num( (nline( (Ptr(2)+3): (Ptr(3)-1) ) ) )/60;         
    if isempty(ImgInfo(i).Geo.Latitude)
       disp('Latitude empty');
       GPSWorks = false;
    end
    if strcmp( nline(Ptr(3)+1), 'S')
        ImgInfo(i).Geo.Latitude = - ImgInfo(i).Geo.Latitude;
    end
      
    % Longitude
    ImgInfo(i).Geo.Longitude = str2num(nline( (Ptr(4)+1):(Ptr(4)+3) ))+str2num(nline( (Ptr(4)+4):(Ptr(5)-1)))/60;       
    if isempty(ImgInfo(i).Geo.Longitude)
       disp('Longtitude empty');
       TargetTime
       PriorTime
       i
       ImgInfo(i).ExifInfo.name
       GPSWorks = false;           
    end
    if strcmp( nline( Ptr(5)+1 ), 'W')
       ImgInfo(i).Geo.Longitude = - ImgInfo(i).Geo.Longitude;
    end
        
    % altitude in meters
    ImgInfo(i).Geo.altitude = str2num(nline( (Ptr(9)+1):(Ptr(10)-1) ));
    if isempty(ImgInfo(i).Geo.altitude)
       disp('altitude empty');
       GPSWorks = false;
    end
      
%     % transform into X Y Z coordinate
%     [X] = Polar2Cartesian(defaultPara, ImgInfo(i).Geo.Longitude, ImgInfo(i).Geo.Latitude, ImgInfo(i).Geo.altitude, true);
%     ImgInfo(i).X_world = X;
    
    % move on on other image
    SumAltitude = SumAltitude + ImgInfo(i).Geo.altitude;
    i = i + 1;
end
SumAltitude = SumAltitude/length(ImgInfo);

for i = 1:length(ImgInfo)
    % Force the same Altitude
    ImgInfo(i).Geo.altitude = SumAltitude;    
end    

return;
