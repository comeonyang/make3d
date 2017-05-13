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
function []=BuildVrmlMetaModel(NewModelFlag, Fdir, Path, InLinePath, R, T, Scale);

if NewModelFlag
	% open the Path file nad write in the first model inline code
	fp = fopen( Path, 'w');
	fprintf(fp, '#VRML V2.0 utf8\n');
	fprintf(fp, 'Transform {\n');
	fprintf(fp, '    translation %f %f %f\n', T);
	fprintf(fp, '    scale %f %f %f\n', Scale);
 	fprintf(fp, '    children [\n');   
    fprintf(fp, '       Transform {\n');
    fprintf(fp, '          rotation %f %f %f %f\n', R);
	fprintf(fp, '          children [\n');
	fprintf(fp, '             Inline {\n');
	fprintf(fp, '                url  "%s"\n',InLinePath);
	fprintf(fp, '             }\n');
	fprintf(fp, '          ]\n');
	fprintf(fp, '       }\n');
	fprintf(fp, '    ]\n');
	fprintf(fp, '}\n');    
	fclose(fp);
else
	% first read to decide to modify or attach new inline code
	if true % Debug: no delete VRML==================
		fp = fopen( Path, 'r')
		out = 1;
		line = 1;
		while out ~=-1
			out = fgets(fp);
			if ~isempty( findstr(out,InLinePath))
				% erase the inline;
		['sed ' num2str(line-6) ',' num2str(line+3) 'd ' Path '> ' Fdir '/temp.wrl']
				system(['sed ' num2str(line-8) ',' num2str(line+5) 'd ' Path '> ' Fdir '/temp.wrl']);
				system(['cp ' Fdir '/temp.wrl ' Path]);
				disp('delete old inline');
				line
				break;	
			end
			line = line +1;
		end
		fclose(fp);
	end % ============================================

	% Start add new inline code
	fp = fopen( Path, 'a+');
	fprintf(fp, 'Transform {\n');
	fprintf(fp, '    translation %f %f %f\n', T);
	fprintf(fp, '    scale %f %f %f\n', Scale);
 	fprintf(fp, '    children [\n');   
	fprintf(fp, '       Transform {\n');
	fprintf(fp, '          rotation %f %f %f %f\n', R);
	fprintf(fp, '          children [\n');
	fprintf(fp, '             Inline {\n');
	fprintf(fp, '                url  "%s"\n',InLinePath);
	fprintf(fp, '             }\n');
	fprintf(fp, '          ]\n');
	fprintf(fp, '       }\n');
	fprintf(fp, '    ]\n');
	fprintf(fp, '}\n');    
	fclose(fp);    
end    
