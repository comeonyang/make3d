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
function [pixel, camCoord, isCorrect] = getPixelCoords(Wcoord,panAngR,tiltAngR, ...
            rigHeight,distCoeffs,camIntMat,UNDISTORT,camExtMat1)
  
isCorrect = true;
camExtMat2 = getCamExtMat(panAngR,tiltAngR,rigHeight); %% its a 4x3 matrix\
% pCamt = camExtMat2*[Wcoord 1]';
% camPanAxisOffset = 0.20;
% camTiltAxisOffset = 0.155;
% pCam(1)=-pCamt(2);
% pCam(2)=-(pCamt(3)-camTiltAxisOffset);
% pCam(3)=pCamt(1)-camPanAxisOffset;
% pCam(4)=pCamt(4);
% 
% pCam=[-pCam(2) pCam(1) pCam(3) 1]'
pCam=camExtMat1*camExtMat2*[Wcoord 1]';
%% check if point behind image plane
if (pCam(3) < 0.0)
    problem_z = 1 %return false;    
    isCorrect = false;
end
camCoord = pCam(1:3);

k1 = -distCoeffs(1);
k2 = -distCoeffs(2);
p1 = -distCoeffs(3);
p2 = -distCoeffs(4);
x = pCam(1)/pCam(3);
y = pCam(2)/pCam(3);
x2 = x*x;
y2 = y*y;
r2 = x2 + y2;
r4 = r2*r2;
xp = x*(1+ k1*r2 + k2*r4) + 2*p1*x*y + p2*(r2+2*x2);
yp = y*(1+ k1*r2 + k2*r4) + 2*p2*x*y + p1*(r2+2*y2);
if(UNDISTORT)
	pixel(1) = camIntMat(1,1)*xp + camIntMat(1,3);
	pixel(2) = camIntMat(2,2)*yp + camIntMat(2,3);
else
	pixel(1) = camIntMat(1,1)*x + camIntMat(1,3);
	pixel(2) = camIntMat(2,2)*y + camIntMat(2,3);
end
if ( (pixel(1) <= 0.0) || (pixel(1) >= 2272) )
    problem_len = 1;    
    isCorrect = false;
end
if ( (pixel(2) <= 0.0) || (pixel(2) >= 1704) )
    problem_wid = 1;    
    isCorrect = false;
end
