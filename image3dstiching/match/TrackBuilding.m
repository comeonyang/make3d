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
function [ track]=TrackBuilding(Matches)

% This function given the Matches between all pairs of images
% generate the track without any inconsistent (delete the track with any inconsistency):
% in that one image observes multiple features in the track

% add image in order from img1, img2, to imgxxx
% When add new image look for matches between added images
% and maintaining the matches

NumImages = size(Matches,1);
track = sparse(0,NumImages);
InConsistMatchesBin = cell(1,NumImages);
NumTrack = size( track,1);
for i = 1:NumImages
	if i ==1 % can build track with only one image
		continue;
	end
	for j = 1:i-1
%		i
%        j
		% Skip if Matches(j,i).Index is empty
		if isempty(Matches(j,i).Index)
			continue;
		end

		% Clean the Matches using InConsistMatchesBin
		NumInconsistMatchesOfJ = size(InConsistMatchesBin{j},1);
		if NumInconsistMatchesOfJ ~=0
			NumMatces = size(Matches(j,i).Index, 2);
			DeleteMark = sum(repmat( Matches(j,i).Index, NumInconsistMatchesOfJ, 1) == repmat(InConsistMatchesBin{j}, 1, NumMatces), 1) ~= 0;
			Matches(j,i).Index(DeleteMark) = [];
			Matches(i,j).Index(DeleteMark) = [];
		end

		NumMatces = size(Matches(j,i).Index, 2);
		% editing existing track
		if NumTrack ~=0
			Ptr = (repmat(track(:,j), 1, NumMatces) == ...
				repmat( Matches(j,i).Index, NumTrack, 1));
			Target = repmat( Matches(i,j).Index, NumTrack, 1);
            NewTarget = sparse(NumTrack, 1);
			NewTarget(sum(Ptr,2)~=0) = sum( Target(Ptr), 2);
            Target = NewTarget;
			Matches_for_new_track_I = setdiff(Matches(i,j).Index, Target');
			% check consistency
			InConsistMark = (Target ~= track(:,i)) & (track(:,i) ~= 0);
			track(:,i) = Target;
			% maintina the InConsistMatchesBin
            if ~all(InConsistMark==0)
                InConsistMatchesBin{j} = union(InConsistMatchesBin{j}, track(InConsistMark,j));
                InConsistMatchesBin{i} = track(InConsistMark,i);
                track(InConsistMark,:) = [];
            end
		else
			Matches_for_new_track_I = Matches(i,j).Index;
		end

		% add new track
		NumNewTrack = length(Matches_for_new_track_I);
		Matches_for_new_track_mark_I = sum( repmat( Matches(i,j).Index, NumNewTrack, 1) ==...
							repmat( Matches_for_new_track_I', 1, NumMatces), 1) ~=0;
		NewTrack = sparse(NumNewTrack, NumImages);
		NewTrack(:,j) = Matches(j,i).Index(Matches_for_new_track_mark_I)';
		NewTrack(:,i) = Matches(i,j).Index(Matches_for_new_track_mark_I)';
		track = [track; NewTrack];
		
		% maintain Num of track
		NumTrack = size( track,1);
	end
end
