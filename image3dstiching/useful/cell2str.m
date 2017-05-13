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
%CELL2STR Convert cell array into evaluable string.
%   B = CELL2STR(C) returns a B such that C = EVAL(B), under the
%   following contraits:
%   - C is composed of numeric arrays or strings.
%   - All of the elements of C are of the same type.
%   - C is a row vector, that is, SIZE(C,1) == 1 and NDIMS(C) = 2.
%
%   See also MAT2STR

% (c)2000 by Cris Luengo

function str = cell2str(c)
disp('hi')

if ~iscell(c)

   if ischar(c)
      str = ['''',c,''''];
   elseif isnumeric(c)
      str = mat2str(c);
   else
      error('Illegal array in input.')
   end

else

   N = length(c);
   if N > 0
      if ischar(c{1})
         str = [c{1}];
         for ii=2:N
            if ~ischar(c{ii})
               error('Inconsistent cell array');
            end
            str = [str c{ii}];
         end
         str = [str];
      elseif isnumeric(c{1})
         str = [mat2str(c{1})];
         for ii=2:N
            if ~isnumeric(c{ii})
               error('Inconsistent cell array');
            end
            str = [str,mat2str(c{ii})];
         end
         str = [str];
      end
   else
      str = '';
   end

end
