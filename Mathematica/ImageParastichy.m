(* ::Package:: *)

makePackage = False;
If[makePackage,
BeginPackage["ImageParastichy`"];
locateParastichyOptions::usage = "Global options for heuristics";
makeTidyMesh::usage = "Delete duplicate points and remove over-large cells";
meshExtendParastichy::usage = "Extend a line without large kinks"
];

If[makePackage,Begin["`Private`"]];



locateParastichyOptions = <|
"DuplicateNearness" ->5.1,
"LargePolygonLength" ->42,
"InitialPeel" -> 0,
"StraightnessWhenAdjacent" ->40,
"StraightnessWhenUnadjacent"-> 30,
"InitialParastichyLength"->  10,
"ParastichyExtensionLength"->  10, (* each way *)
"MinimumThreadLength"-> 3,
"FamilyGrowthSize" -> 20
|>; 
	
makeTidyMesh[seedCentres_] :=
	Module[{dMeshRaw, lineLength, lengths, pointsToDrop, x, meshxy, toDrop, dMesh, res, boundary},
		dMeshRaw = DelaunayMesh[seedCentres];
	
	lengths = AssociationThread[MeshCells[dMeshRaw, 1], Map[meshLineLength, MeshPrimitives[dMeshRaw, 1]]];
	(* find the short line pairs and take the first cell index *) 
	pointsToDrop =
	(
	Keys @ Select[lengths, # < locateParastichyOptions["DuplicateNearness"]&] /. Line[{x_, y_}] -> x
	);
	toDrop = Map[First @ Position[seedCentres, meshCoordinates[dMeshRaw,#]]&, pointsToDrop];
	(* should be in point order anyway *)
	Print["Dropping ", Length[toDrop], " points as duplicates"];
	dMesh = DelaunayMesh[Delete[seedCentres, toDrop]];
	dMesh = meshDeletePolygons[dMesh, meshLongCells[dMesh, locateParastichyOptions["LargePolygonLength"]]];
	dMesh = meshPeel[dMesh, locateParastichyOptions["InitialPeel"]];
	boundary = SortBy[ConnectedMeshComponents @ RegionBoundary[dMesh], Length[MeshCells[#, 0]&]];
	res = <|
			"Mesh" -> dMesh
			, "InnerBoundary" -> meshCellsOfMesh1InMesh2[boundary[[1]], dMesh]
			, "OuterBoundary" -> meshCellsOfMesh1InMesh2[boundary[[2]], dMesh]
			, "Adjacency" -> meshAdjacencyAssociation[dMesh]
		|>;
		res
	];	



(* private functions *)


meshLineLength[Line[{p1_, p2_}]] := Norm[p1 - p2];
meshCoordinates[mesh_MeshRegion,ix_] := Block[{x},MeshPrimitives[mesh, {0, ix}] /. Point[x_] -> x];
meshCoordinates[meshAssociation_Association,ix_] := meshCoordinates[meshAssociation["Mesh"],ix];
meshCoordinates[mesh_] := Block[{x},MeshPrimitives[mesh, 0] /. Point[x_] -> x];

meshDeletePolygons[mesh_,deletedpolys_] :=   
	MeshRegion[meshCoordinates[mesh],Complement[ MeshCells[mesh,2],deletedpolys]];

meshLinesWithLength[mesh_MeshRegion] := AssociationThread[MeshCells[mesh,1],Map[meshLineLength,MeshPrimitives[mesh,1]]];

(* for tidying the initial mesh *)

meshLongCells::usage = "Find cells in which one of the sides has length> length ";
meshLongCells[mesh_,length_] := Module[{polys,x,lengths,polyMembers,longLinePairs,polyHasLinePair,polyHasAnyPairs},
	longLinePairs = Keys@Select[meshLinesWithLength[mesh], # > length&] /. Line[x_]->x;
	polys = MeshCells[mesh,2];
	polyMembers  =AssociationThread[polys,polys/. Polygon[x_]->x];
	polyHasLinePair[Polygon[x_],pair_] := Length[Intersection[x,pair]]>= 2;
	polyHasAnyPairs[p_] :=Or@@Map[polyHasLinePair[p,#]&,longLinePairs];
	Keys@Select[Association@Map[#->polyHasAnyPairs[#]&,polys],#&]
];

(* peel the outer boundary of cells off n times *)
meshPeel[mesh_,n_] := Module[{doapeel},
	doapeel[meshp_] := MeshRegion[MeshPrimitives[meshp,0]/. Point[x_]->x,meshPolygonsNotFormingBoundary[meshp]];
	Nest[doapeel,mesh,n]
];
meshPolygonsNotFormingBoundary[mesh_] := Module[{boundaryNodes,meshPolygons,polyFormsBoundary},
	boundaryNodes = boundaryPointIndexes[mesh];
	meshPolygons = MeshCells[mesh,2];
	polyFormsBoundary[Polygon[threeNodeList_]] := Length[Intersection[threeNodeList,boundaryNodes]]>= 2;
	Select[meshPolygons,Not @* polyFormsBoundary]
];

boundaryPointIndexes[mesh_] := Module[{boundaryNodesXY,meshXY},
(* replace by cellsOfMesh1InMesh2[RegionBoundary[mesh], mesh] just for clarity *)
	boundaryNodesXY = MeshPrimitives[RegionBoundary[mesh],0];
	meshXY = MeshPrimitives[mesh,0];
	Map[First@FirstPosition[meshXY,#]&,boundaryNodesXY]
];

meshCellsOfMesh1InMesh2[mesh1_, mesh2_] :=
	Module[{pointxy},
		pointxy = MeshPrimitives[mesh1, 0];
		Flatten @ Map[FirstPosition[MeshPrimitives[mesh2, 0], #]&, pointxy]
	];

meshAdjacencyAssociation[mesh_] := Module[{paraPossibles},
	paraPossibles =  Flatten[Apply[List,MeshCells[mesh,1],{1}],1];
	paraPossibles = DeleteDuplicates@Join[paraPossibles,Map[Reverse,paraPossibles]];
	GroupBy[paraPossibles,First,Map[Last,#]&]
];

If[makePackage,End[]]; (* Private *)

If[makePackage,EndPackage[]];





(* ::Input::Initialization:: *)

createParastichyFamily[meshAssociation_,starter_,family_:1] := Module[{nextpara,paraList,mpf,pstart},
pstart = parastichyStarter[meshAssociation,starter];If[Length[pstart]<2,Return[Missing["Can't make starter from ", starter]]];
mpf = {};
mpf = addToParastichyFamily[mpf,pstart,family];
mpf = Nest[findAdjacentThreads[meshAssociation,#]&,mpf,
locateParastichyOptions["FamilyGrowthSize" ]];
mpf = tidyParastichyFamily[meshAssociation,mpf];
mpf = renumberParastichyFamily[meshAssociation,mpf];
mpf
];

tidyParastichyFamily[meshAssociation_,parastichyFamily_] := Module[{family,ix,jx,res,overlap,newParastichy,lastIndex,firstTail,secondHead},
res= parastichyFamily;family = First[parastichyFamily]["Family"];
For[ix=1,ix<Length[parastichyFamily],ix++,
For[jx = ix+1, jx<=Length[parastichyFamily],jx++,
overlap  = Intersection[parastichyFamily[[ix]]["Members"],parastichyFamily[[jx]]["Members"]];
If[Length[overlap]>0,
res = mergeParastichySibs[meshAssociation,parastichyFamily,res,{ix,jx}];
If[MissingQ[res],Echo[{ix,jx,res, "overlapping failed"},"tPF"];Return[parastichyFamily]];
];
firstTail = Last[parastichyFamily[[ix]]["Members"]];
secondHead =  First[parastichyFamily[[jx]]["Members"]];
If[MemberQ[meshAssociation["Adjacency"][firstTail],secondHead],
res = mergeParastichySibs[meshAssociation,parastichyFamily,res,{ix,jx}];
If[MissingQ[res],Echo[{{ix,jx},{parastichyFamily[[ix]],parastichyFamily[[jx]]},firstTail,secondHead, "joining failed"},"tPF"];Return[parastichyFamily]];]
]];
res
];



renumberParastichyFamily[meshAssociation_,parastichyFamily_] := Module[{g,path},
paraAdjacency[parastichyFamilyMember_] := DeleteDuplicates[Union@@Map[meshAssociation["Adjacency"],parastichyFamilyMember["Members"]]];
paraAdjacencyCount[ix_,ix_] := 0;
paraAdjacencyCount[ix_,jx_] := Length@Intersection[
paraAdjacency[parastichyFamily[[ix]]],
parastichyFamily[[jx]]["Members"]
];
paraAdjacencies[ix_] := Keys@Take[Reverse@Sort@Association@Table[ 
UndirectedEdge[parastichyFamily[[ix]]["Index"],parastichyFamily[[jx]]["Index"]]-> paraAdjacencyCount[ix,jx],{jx,Length[parastichyFamily]}],2];
atable = Flatten@Table[paraAdjacencies[ix],{ix,Length[parastichyFamily]}];
g = Graph[Map[#["Index"]&,parastichyFamily],atable,VertexLabels->"Name"];
path = findSpanningPathEitherWay[g];
res = {};family = First[parastichyFamily]["Family"];family=ToString[family]~~"A";
For[i=1,i<= Length[path],i++,
para = Query[SelectFirst[#["Index"]==path[[i]]&]]@parastichyFamily;
res = addToParastichyFamily[res,para["Members"],family];
];
res
];

mergeParastichySibs[meshAssociation_,masterFamily_,newFamily_,{ix_,jx_}] := Module[{family,res,newParastichy
,ixxi,jxxj,decho},
decho[x_] := If[{ix,jx}=={5,-17},Echo[x,"mPS"],x];
res = newFamily;
family = First[masterFamily]["Family"];
newParastichy = Union[masterFamily[[ix]]["Members"],masterFamily[[jx]]["Members"]];
newParastichy = decho@makeDirectedParastichy[meshAssociation,newParastichy];

ixxi  = FirstPosition[res,SelectFirst[res, #["Index"]==masterFamily[[ix]]["Index"]&]];
If[!MissingQ[ixxi],res = Drop[res,ixxi]];
jxxj =  FirstPosition[res,SelectFirst[res, #["Index"]==masterFamily[[jx]]["Index"]&]];
If[!MissingQ[jxxj],res = Drop[res,jxxj]];
res = addToParastichyFamily[res,newParastichy,family];
res
];


addToParastichyFamily[mpf_,parastichy_,family_] := Module[{ix,res,overlaps,maxIndex},
res = mpf;
maxIndex = If[Length[mpf]==0,0,Max[Map[#["Index"]&,mpf]]];
res = Append[res,
<| "Family"->family,"Index"-> maxIndex+1,"Members"->parastichy , "Head"-> First@parastichy|>];
res
];

parastichyStarter[meshAssociation_, starter_] :=
  Module[{inner, parastichyPoints},
    inner = meshExtendParastichy[meshAssociation, starter, locateParastichyOptions["InitialParastichyLength"], {}, {}];
   inner
  ];
 
meshExtendParastichy[meshAssociation_,starter_,length_,avoidPoints_,adjacentParastichy_] :=Module[{i,res,np,decho},
decho[x_] := Echo[x,"mEP"];
res = starter;

For[ i=1,i<= length,i++,
np = nextParastichyPoint[meshAssociation,res,avoidPoints,adjacentParastichy];
If[MissingQ[np],Break[]];res= Append[res,np];If[MemberQ[meshAssociation["OuterBoundary"],np],Break[]];If[MemberQ[meshAssociation["InnerBoundary"],np],Break[]];
];
res
];


nextParastichyPoint[meshAssociation_,parastichy_,avoidPoints_,adjacentParastichy_] := Module[{candidateContinues,straightnessAngle,nextStraightest,debugTest,adjacentParastichyEnds,candidateInfo,adjacentParastichyCentre,
adjacentCandidates,res,decho},
debugTest = MemberQ[{2,10},-Last@parastichy];
decho[x_] := If[debugTest,Echo[x,"nPP"],x];

If[Length[parastichy]==0,Return[Missing["No parastichy to extend"]]];If[MissingQ[parastichy]==0,Return[Missing["Missing parastichy"]]];

straightnessAngle = If[Length[adjacentParastichy]>0,
locateParastichyOptions["StraightnessWhenAdjacent"],
locateParastichyOptions["StraightnessWhenUnadjacent"]
];
adjacentParastichyEnds  = If[Length[adjacentParastichy]==0,{},Union@@Map[meshAssociation["Adjacency"][#]&, {First[adjacentParastichy],Last[adjacentParastichy]}]];adjacentParastichyCentre   = Complement[adjacentParastichy,adjacentParastichyEnds];

decho[parastichy];
candidateContinues = meshAssociation["Adjacency"][Last@parastichy];

candidateContinues = Complement[candidateContinues,avoidPoints];candidateInfo = KeyValueMap[ <| "Node"->#1,"Deviation"-> #2|> &]@orderByStraightness[meshAssociation,parastichy,candidateContinues];candidateInfo = Map[
Append[#, "Adjacency"-> meshAssociation["Adjacency"][#["Node"]]]&,candidateInfo];candidateInfo = Map[Append[#, "NextToAdjacentCentre"->  
Length[Intersection[#["Adjacency"],adjacentParastichyCentre]]>0]&,candidateInfo];candidateInfo = Map[Append[#, "NextToAdjacentEnds"->  
Length[Intersection[#["Adjacency"],adjacentParastichyEnds]]>0]&,candidateInfo];candidateInfo = Map[
Append[#, "NextToAdjacent"->  
Length[Intersection[#["Adjacency"],adjacentParastichy]]>0]&,candidateInfo];candidateInfo = Map[
Append[#, "OnAdjacent"-> 
 Length[Intersection[{#["Node"]},adjacentParastichy]]>0]&,candidateInfo];
candidateInfo = Map[
Append[#, "Straightish"-> isStraightish[#]]&,candidateInfo];
(* screen out any with too much deviation anyway *) 

candidateInfo =  Query[Select[  #["Straightish"]& ]]@ candidateInfo;If[Length[candidateInfo]==0,
decho@"End of the line" ;Return[Missing["End of the line"]]];If[First[candidateInfo]["OnAdjacent"],decho@{"terminating at ",candidateInfo, "because of collision with adjacent"};Return[Missing["Collision"]]
];

(* do we have any that preserve adjacency *) adjacentCandidates =  Query[Select[  #["NextToAdjacent"]& ]]@candidateInfo;If[debugTest,Print[adjacentCandidates ]];If[Length[adjacentCandidates]>0,
res = First[adjacentCandidates]["Node"];
If[debugTest,Print[res]];
Return[res]];

(* otherwise we take a nonadjacent (which may have been straighter *)res = First[candidateInfo]["Node"];If[debugTest,Print[res]];Return[res];
];

isStraightish[cInfo_] := (
Abs[cInfo["Deviation"]]< 
If[cInfo["NextToAdjacentCentre"] || cInfo["NextToAdjacentEnds"], 
locateParastichyOptions["StraightnessWhenAdjacent"],
locateParastichyOptions["StraightnessWhenUnadjacent"]
]
);

orderByStraightness[meshAssociation_,parastichy_,candidateContinues_] := Module[{paraPossiblePairs,paraPossibleAngles,lastAngle,paraAngles,deviationsNext},paraPossiblePairs = Map[{parastichy[[-1]],#}&,candidateContinues];paraPossibleAngles = Association@Map[#[[2]]->meshLineAngle[meshAssociation,#]&,paraPossiblePairs];lastAngle = meshLineAngle[meshAssociation,{parastichy[[-2]],parastichy[[-1]]}];paraAngles = Map[anglePrincipal[#-lastAngle]&,paraPossibleAngles];deviationsNext = Sort@Abs[paraAngles];
deviationsNext
];

meshLineDeviation[meshAssociation_,{ix1_,ix2_,ix3_}] := Module[{p1p2,p2p3,first,second,pAngle,res},
(* in [-180,180 *)p1p2 = MeshPrimitives[meshAssociation["Mesh"],{0,{ix1,ix2}}];pAngle[{Point[{x1_,y1_}],Point[{x2_,y2_}]}] := (360/(2\[Pi])) ArcTan[x2-x1,y2-y1];first = pAngle[p1p2];
p2p3 = MeshPrimitives[meshAssociation["Mesh"],{0,{ix2,ix3}}];
second = pAngle[p2p3];
anglePrincipal[first-second]
];

meshLineAngle[meshAssociation_,{ix1_,ix2_}] := Module[{p1p2,pAngle,res},
(* in [-180,180 *)
p1p2 = MeshPrimitives[meshAssociation["Mesh"],{0,{ix1,ix2}}];
pAngle[{Point[{x1_,y1_}],Point[{x2_,y2_}]}] := (360/(2\[Pi])) ArcTan[x2-x1,y2-y1];
pAngle[p1p2]
];
anglePrincipal[angle_] :=  angle - 360 Round[angle/360]; (* in -180 < angle < 180 *) 



(* ::Input::Initialization:: *)
findAdjacentThreads[meshAssociation_,parastichyFamily_] := Module[{res,atp,family,ix,avoidPoints,extendThread,decho},
atp = adjacentThreadsToParastichyFamily[meshAssociation,parastichyFamily];
decho[x_] := If[debugTest,Echo[x,"fAT"],x];
debugTest = MemberQ[{},Last[parastichyFamily]["Index"]];
decho[Last[parastichyFamily]["Index"]];
decho[atp];
(* a list of directed paths *)
family = Last[parastichyFamily]["Family"];
res = parastichyFamily;
For[ix=1,ix<= Length[atp],ix++,
avoidPoints = nodesInParastichyFamily[parastichyFamily];
extendThread = meshExtendParastichy[meshAssociation,atp[[ix]]
,locateParastichyOptions["ParastichyExtensionLength"], avoidPoints,{}];
avoidPoints = Union[avoidPoints,extendThread];
extendThread = meshExtendParastichy[meshAssociation,Reverse@extendThread
,locateParastichyOptions["ParastichyExtensionLength"], avoidPoints,{}];
extendThread  = makeDirectedParastichy[meshAssociation,extendThread];
res = addToParastichyFamily[res,extendThread,family] ; 
];
res
];

centralMembers[parastichyFamilyMember_] := Module[{members},
members = parastichyFamilyMember["Members"];
If[Length[members]==0,Return[Missing["Error: parastichy with no members"]]];
If[Length[members]<3,Return[{}]];
Drop[Drop[members,-1],1]
];

nodesInParastichyFamily[parastichyFamily_] :=  DeleteDuplicates[Union@@Map[#["Members"]&,parastichyFamily]];

nonEndNodesInParastichyFamily[parastichyFamily_] :=  DeleteDuplicates[Union@@Map[centralMembers,parastichyFamily]];

adjacentNodesToParastichyFamily[meshAssociation_,parastichyFamily_]:=  Module[{familyPoints,allAdjacentPoints},
familyPoints  =  nonEndNodesInParastichyFamily@parastichyFamily;
allAdjacentPoints = Union @@ Map[meshAssociation["Adjacency"][#]&, familyPoints];
allAdjacentPoints = DeleteDuplicates @allAdjacentPoints ;
allAdjacentPoints= Complement[allAdjacentPoints, nodesInParastichyFamily[parastichyFamily] ];
allAdjacentPoints
];

adjacentThreadsToParastichyFamily[meshAssociation_,parastichyFamily_]:=  Module[{allAdjacentPoints,allAdjacentGraph,allAdjacentComponents,adjacentThreads},
allAdjacentPoints =adjacentNodesToParastichyFamily[meshAssociation,parastichyFamily];
If[Length[allAdjacentPoints] < 2,
			Return[{}]
];
allAdjacentGraph= adjacencyGraph[meshAssociation, allAdjacentPoints];
allAdjacentComponents =  ConnectedGraphComponents[allAdjacentGraph];
allAdjacentComponents  =Select[allAdjacentComponents, VertexCount[#] >= locateParastichyOptions["MinimumThreadLength"]&];

allAdjacentComponents  =Map[makeDirectedPath[meshAssociation,#]&,allAdjacentComponents];

allAdjacentComponents    = allAdjacentComponents /. Missing[_]->Nothing[];
allAdjacentComponents  = Map[splitAtKinks[meshAssociation,#]&,allAdjacentComponents];
allAdjacentComponents;
allAdjacentComponents  = Flatten @allAdjacentComponents;
allAdjacentComponents    = allAdjacentComponents /. Missing[_]->Nothing[];
adjacentThreads  =Map[TopologicalSort,allAdjacentComponents];
adjacentThreads
];

vertexCount[graph_,vertex_] := Length[IncidenceList[graph,vertex]];

removeCyclicPoints[component_] :=Module[{stree,incidence,edgeCount,branchPoints,streeComponents,res},
stree = FindSpanningTree[component];
incidence  =  AssociationThread[VertexList[stree],
Map[IncidenceList[stree,#]&,VertexList[stree]]
];
edgeCount  =  Map[Length,incidence];
branchPoints = Select[edgeCount,#>2&];
stree =VertexDelete[component,Keys@branchPoints];
streeComponents = ConnectedGraphComponents[stree];
streeComponents = SortBy[streeComponents,VertexCount];
res = Last@streeComponents;
res

];

findSpanningPathEitherWay[undirectedGraph_] := Module[{acyclic,endPoints,path,decho},
decho[x_] := Echo[x,"fSPEW"];
acyclic = FindSpanningTree[undirectedGraph];

endPoints = Select[VertexList[acyclic],vertexCount[acyclic,#]==1&];
If[Length[endPoints]<1,Return[Missing["Missing path ends"]]];

If[Length[endPoints]==1,
decho[{"No spanning tree for ",acyclic}];
Return[Missing["Couldn't construct spanning tree"]]
];

If[Length[endPoints]>2,endPoints=Take[endPoints,2]];

path = FindPath[acyclic,First[endPoints],Last[endPoints]];
If[Length[path]==0,
path = FindPath[acyclic,Last[endPoints],First[endPoints]]
];

If[Length[path]==0,
Return[Missing["No spanning path"]]];

path = First@path;
path
];

findSpanningPath[meshAssociation_,undirectedGraph_] := Module[{acyclic,endPoints,path,decho},
decho[x_] := Echo[x,"fSPd"];
acyclic = removeCyclicPoints[undirectedGraph];

endPoints = Select[VertexList[acyclic],vertexCount[acyclic,#]==1&];
endPoints = SortBy[endPoints,-meshPointRadius[meshAssociation,#]&];

If[Length[endPoints]<1,Return[Missing["Missing path ends"]]];

If[Length[endPoints]==1,
decho[{"No spanning tree for ",acyclic}];
Return[Missing["Couldn't construct spanning tree"]]];

If[Length[endPoints]>2,endPoints=Take[endPoints,2]];

path = FindPath[acyclic,First[endPoints],Last[endPoints]];
If[Length[path]==0,
path = FindPath[acyclic,Last[endPoints],First[endPoints]]
];

If[Length[path]==0,
Return[Missing["No spanning path"]]];
Clear[decho];
path = First@path;
path
];

makeDirectedParastichy[meshAssociation_,parastichy_] := Module[{edgeList,path},
edgeList = Map[{List,{#},Intersection[parastichy,meshAssociation["Adjacency"][#]]}&,parastichy];
edgeList = DeleteDuplicates[Sort /@Flatten[Map[Apply[Outer,#]&,edgeList],2]];
edgeList = Map[ #[[1]] \[UndirectedEdge] #[[2]] &, edgeList];
path  = findSpanningPath[meshAssociation, Graph[edgeList]];
path
];

makeDirectedPath[meshAssociation_,component_] :=Module[{internalVertices,endPoints,path,acyclic},
(* given an undirected graph which is mainly a single path , return the directed graph with outer point at head *)
acyclic = removeCyclicPoints[component];
path = findSpanningPath[meshAssociation,acyclic];
If[MissingQ[path],Return[path]];
path = Partition[path,2,1];
res= EdgeDelete[acyclic,EdgeList[acyclic]];
edgesToAdd = DirectedEdge @@@ path;
res = EdgeAdd[res,DirectedEdge @@@ path];
res
];
directedPathHead[component_] :=  First@TopologicalSort[component];


meshPointRadius[mesh_MeshRegion,ix_] := Module[{c,p},
c = RegionCentroid[mesh];
p = meshCoordinates[mesh,ix];
Norm[(c-p)]
];

meshPointRadius[meshAssociation_Association,ix_] := meshPointRadius[meshAssociation["Mesh"],ix]

vertexAngle[meshAssociation_,component_,vertex_] :=Module[{incidence},
incidence = IncidenceList[component,vertex];
If[Length[incidence]==0,Return[Missing["Can't find incidence"]]];
If[Length[incidence]==1,
Return[0]];
If[Length[incidence] > 2,Return[Missing["Too many incidences"]]];
inPair = First@Cases[incidence, _ \[DirectedEdge] vertex];
outPair = First@Cases[incidence, vertex \[DirectedEdge] _];
ix123 = {First[inPair],First[outPair],Last[outPair]};
meshLineDeviation[meshAssociation,ix123]
];

splitAtKinks[meshAssociation_,component_] :=
(* component is an unforked directed path with one element *) 
Module[{angles},

angles = Association@Map[#->vertexAngle[meshAssociation,component,#]&,
VertexList[component]];
kinks = Keys@Select[angles,Abs[#] >
 locateParastichyOptions["StraightnessWhenAdjacent"]
 &];
edgesToDrop = Map[
Cases[IncidenceList[component,#],#\[DirectedEdge] _]&,kinks];
res = component;
res = EdgeDelete[res,Flatten@edgesToDrop];
res  =connectedDirectedGraphComponents[res];
res = Select[res,VertexCount[#]>  locateParastichyOptions["MinimumThreadLength"]&];

Return[res];

];

edeleteEdge[g_Graph,edge_] :=If[EdgeQ[g,edge],EdgeDelete[g,edge],g];

edelete[comp_,edges_] := Module[{g,i},
g=comp;
For[i=1,i<= Length[edges],i++,
g = edeleteEdge[g,edges[[i]]]
];
g
];

connectedDirectedGraphComponents[g_] := Module[{edges,flip,res},
edges = EdgeList[g];
flip[ a_ \[DirectedEdge] b_ ] := b \[DirectedEdge] a;
edges = Map[flip,edges];
If[MissingQ[edges],Echo[" Missing edges ", Stack[]]];
res = EdgeAdd[g,edges];
Map[edelete[#,edges]&,
ConnectedGraphComponents[res]]
];


adjacencyGraph[meshAssociation_,vertexes_] := Module[{g,edges,edge,adjacencies,vertexCoordinates},

adjacencies= KeyTake[meshAssociation["Adjacency"],vertexes];adjacencies= Map[Intersection[vertexes,#]&,adjacencies];
makeSortedEdge[v1_,v2_] := UndirectedEdge @@ Sort[{v1,v2}];
edge[k_,v_] := Map[ makeSortedEdge[k,#]&,v];edges =DeleteDuplicates@Flatten[KeyValueMap[edge,adjacencies]];

vertexCoordinates := Map[meshCoordinates[meshAssociation["Mesh"],{0,#}]&,vertexes] ;g  = Graph[vertexes,edges,VertexLabels -> "Name",
 VertexCoordinates -> vertexCoordinates];
g
];

