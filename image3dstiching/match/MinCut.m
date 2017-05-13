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
% performs Min Cut
% input - list of vertices to keep togeather and weighted      CLICK !!! <=====
% output -
% MinCutGroupsList - two lists of verices, SECOND one contains the sourve vertives
% MinCutWeight - weight on cut
function [MinCutGroupsList, MinCutWeight] = MinCut(SourceVertices, WeightedGraph)

GraphDim = size(WeightedGraph,1);
SourceVertices = SourceVertices(find(SourceVertices));

% remove self edges and ZEROed ones
for ii = 1:GraphDim
    WeightedGraph(ii,ii) = inf;
end
WeightedGraph(find(WeightedGraph == 0)) = inf;

%Merge all Source Vrtices to one, so they'll be unbreakable, descending order is VITAL!!!
SourceVertices = sort(SourceVertices);
GroupsList = zeros(GraphDim);   %each row are the vertices melted into one vertex in the table.
GroupsList(:,1) = 1:GraphDim;
for ii=length(SourceVertices):-1:2; 
    [WeightedGraph,GroupsList] = MeltTwoVertices(SourceVertices(1),SourceVertices(ii),WeightedGraph,GroupsList);
end
Split = GroupsList(:,1);

[MinCutGroupsList_L, MinCutWeight] = MinCutNoSeed(WeightedGraph);
for ii = 1:2
    MinCutGroupsList(ii,:) = Local2GlobalIndices(MinCutGroupsList_L(ii,:), Split);
end

if (length(find(MinCutGroupsList(1,:) == SourceVertices(1))) == 1)
    SeedLocation = 1;
else
    SeedLocation = 2;    
end
MinCutGroupsList_withSeed = [MinCutGroupsList(SeedLocation,find(MinCutGroupsList(SeedLocation,:))) SourceVertices(2:length(SourceVertices))];
MinCutGroupsList_withSeed = sort(MinCutGroupsList_withSeed);
MinCutGroupsList_withSeed = [MinCutGroupsList_withSeed zeros(1,GraphDim - length(MinCutGroupsList_withSeed))];

MinCutGroupsList_NoSeed = MinCutGroupsList(3 - SeedLocation,find(MinCutGroupsList(3 - SeedLocation,:)));
MinCutGroupsList_NoSeed = sort(MinCutGroupsList_NoSeed);
MinCutGroupsList_NoSeed = [MinCutGroupsList_NoSeed zeros(1,GraphDim - length(MinCutGroupsList_NoSeed))];

MinCutGroupsList = [MinCutGroupsList_NoSeed ; MinCutGroupsList_withSeed];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [MinCutGroupsList, MinCutWeight] = MinCutNoSeed(WeightedGraph)
GraphDim = size(WeightedGraph,1);
GroupsList = zeros(GraphDim);   
GroupsList(:,1) = 1:GraphDim;

MinCutWeight = inf;
MinCutGroup = [];
for ii = 1:GraphDim-1
    [OneBefore, LastVertex] = MinimumCutPhase(WeightedGraph);
    if OneBefore == -1 %Graph is not connected. LastVertex is a group of vertices not connected to Vertex 1
            MinCutGroup_L = LastVertex(find(LastVertex)); clear LastVertex; %it's not the last vertex

            MinCutGroupsList = [];
            for jj = 1:length(MinCutGroup_L);
                MinCutGroup_temp = GroupsList(MinCutGroup_L(jj)); 
                MinCutGroup_temp = MinCutGroup_temp(find(MinCutGroup_temp));
                MinCutGroupsList = [MinCutGroupsList MinCutGroup_temp];
            end
            MinCutGroupsList = [MinCutGroupsList zeros(1,GraphDim - length(MinCutGroupsList))];
                     
            jj = 1;
            for kk=1:GraphDim
                if length(find(MinCutGroupsList(1,:) == kk)) == 0
                    MinCutGroupsList(2 ,jj) = kk;
                    jj = jj + 1;
                end
            end
            MinCutWeight = 0;
            return
    end %of: If graph is not connected
    Edges = WeightedGraph(LastVertex,:);
    Edges = Edges(isfinite(Edges));
    MinCutPhaseWeight = sum(Edges);
    if MinCutPhaseWeight < MinCutWeight
        MinCutWeight = MinCutPhaseWeight;
        MinCutGroup = GroupsList(LastVertex,:);
        MinCutGroup = MinCutGroup(find(MinCutGroup));
    end
    [WeightedGraph,GroupsList] = MeltTwoVertices(LastVertex,OneBefore,WeightedGraph,GroupsList);
end

MinCutGroup = sort(MinCutGroup);
MinCutGroupsList = [MinCutGroup zeros(1,GraphDim - length(MinCutGroup))];

jj = 1;
for kk=1:GraphDim
    if length(find(MinCutGroup(1,:) == kk)) == 0
        MinCutGroupsList(2 ,jj) = kk;
        jj = jj + 1;
    end
end

% This function takes V1 and V2 as vertexes in the given graph and MERGES
% THEM INTO V1 !!
function        [UpdatedGraph,GroupsList] = MeltTwoVertices(V1,V2,WeightedGraph,GroupsList)
t = min(V1,V2);
V2 = max(V1,V2);
V1 = t;

GraphDim = size(WeightedGraph,1);
UpdatedGraph = WeightedGraph;

Mask1 = isinf(WeightedGraph(V1,:) );
Mask2 = isinf(WeightedGraph(V2,:) );
UpdatedGraph(V1,Mask1) = 0;
UpdatedGraph(V2,Mask2) = 0;
infMask = zeros(1,size(Mask1,2));
infMask(find(Mask1 & Mask2)) = inf;
UpdatedGraph(V1,:)  =  UpdatedGraph(V1,:)  + UpdatedGraph(V2,:) + infMask;
UpdatedGraph(:,V1) = UpdatedGraph(V1,:)';
UpdatedGraph = [UpdatedGraph(1:V2-1,:) ; UpdatedGraph(V2+1:GraphDim,:)]; %remove second vertex from graph
UpdatedGraph = [UpdatedGraph(:,1:V2-1)  UpdatedGraph(:,V2+1:GraphDim)];
UpdatedGraph(V1,V1) = inf;                                                                                                              % group-group distance

V1list = GroupsList(V1,find( GroupsList(V1,:) ) );
V2list = GroupsList(V2,find( GroupsList(V2,:) ) );
GroupsList(V1,:) = [V1list V2list zeros(1,size( GroupsList,2)- length(V1list) - length(V2list)) ]; %shorten grouplist
GroupsList = [GroupsList(1:V2-1,:) ;GroupsList(V2+1:GraphDim,:) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return [-1, B  ] in case of Unconnected Graph when B is a subgraph(s)
% that are not connected to Vertex 1
function [OneBefore, LastVertex] = MinimumCutPhase(WeightedGraph)
GraphDim = size(WeightedGraph,1);
GroupsList = zeros(GraphDim);
GroupsList(:,1) = 1:GraphDim;

if size(WeightedGraph,1) > 2
    FarestVertexGroup = [0];
    for ii = 1:GraphDim-1
        OneBefore         = FarestVertexGroup(1);
        PossibleVertices = WeightedGraph(1,1:size(WeightedGraph,2));
        PossibleVertices(isinf(PossibleVertices)) = 0;
        FarestVertex = min(find(PossibleVertices == max(PossibleVertices)));
        if FarestVertex == 1 %In case the Graph is not connected
            OneBefore = -1;
            LastVertex = GroupsList(1,:);
            return
        end
        FarestVertexGroup = GroupsList(FarestVertex,:);
        [WeightedGraph,GroupsList] = MeltTwoVertices(1,FarestVertex,WeightedGraph,GroupsList);
    end
    LastVertex = FarestVertexGroup(1);
else
    OneBefore    = 1;
    LastVertex = 2;
end
