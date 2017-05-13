# *  This code was used in the following articles:
# *  [1] Learning 3-D Scene Structure from a Single Still Image, 
# *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
# *      In ICCV workshop on 3D Representation for Recognition (3dRR-07), 2007.
# *      (best paper)
# *  [2] 3-D Reconstruction from Sparse Views using Monocular Vision, 
# *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
# *      In ICCV workshop on Virtual Representations and Modeling 
# *      of Large-scale environments (VRML), 2007. 
# *  [3] 3-D Depth Reconstruction from a Single Still Image, 
# *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
# *      International Journal of Computer Vision (IJCV), Aug 2007. 
# *  [6] Learning Depth from Single Monocular Images, 
# *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
# *      In Neural Information Processing Systems (NIPS) 18, 2005.
# *
# *  These articles are available at:
# *  http://make3d.stanford.edu/publications
# * 
# *  We request that you cite the papers [1], [3] and [6] in any of
# *  your reports that uses this code. 
# *  Further, if you use the code in image3dstiching/ (multiple image version),
# *  then please cite [2].
# *  
# *  If you use the code in third_party/, then PLEASE CITE and follow the
# *  LICENSE OF THE CORRESPONDING THIRD PARTY CODE.
# *
# *  Finally, this code is for non-commercial use only.  For further 
# *  information and to obtain a copy of the license, see 
# *
# *  http://make3d.stanford.edu/publications/code
# *
# *  Also, the software distributed under the License is distributed on an 
# * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# *  express or implied.   See the License for the specific language governing 
# *  permissions and limitations under the License.
# *
# */
#!/bin/bash

dir=$1
img1=$2
img2=$3
Type=$4
abs_thre=$5
ratio_thre=$6

if [ "$Type" == "_" ]; then
	Type="";
	echo 'Default Type';
fi
#echo "Creating surf matches for images of image1 and image2 in dir"

mkdir -p $dir/surf_matches
mkdir -p $dir/match_images

if [ "$img1" != "$img2" ]; then
#   echo "Matching $img1 to $img2"

                    #gen="-s" #exclude non-symetric matches,
   gen="-s -im1 $dir/pgm/$img1.pgm -im2 $dir/pgm/$img2.pgm \
        -o $dir/match_images/"$img1"_"$img2".pgm -c"
#        -o $dir/match_images/$img1-$img2.pgm -c"
                    #gen="-im1 $dir/pgm/$img1.pgm -im2 $dir/pgm/$img2.pgm \
                    #  -o $dir/match_images/$img1-$img2.pgm -c"

                #../../third_party/SURF-V1.0.9/match.ln -k1 $dir/surf/$img1.surf \
#   ../../third_party/SURF-V1.0.8/symmatch.ln -k1 $dir/surf/$img1.surf_64 \
#        -k2 $dir/surf/$img2.surf_64 $gen | sed -e 's/ Matched feature //g' | \
   ../../third_party/SURF-V1.0.8/symmatch_thre.ln \
	-abs $abs_thre -ratio $ratio_thre \
	-k1 $dir/surf/$img1.surf$Type \
        -k2 $dir/surf/$img2.surf$Type \
	$gen | sed -e 's/ Matched feature //g' | \
        sed -e 's/ in image 1 with feature / /g' | \
        sed -e 's/ in image 2\.//g' | sed '$d' > $dir/surf_matches/$img1-$img2.match$Type"_"$abs_thre"_"$ratio_thre
#        sed -e 's/ in image 2\.//g' | sed '$d' > $dir/surf_matches/$img1-$img2.match$Type
fi

