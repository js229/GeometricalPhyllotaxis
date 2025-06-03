(* ::Package:: *)

BeginPackage["DiskStacking`"];


executeRun::usage = "Principal run code";


runFromLattice::usage = "Create run initial condition from a Lattice object";
restartRunFromRun::usage = "Restart from previous run";

pruneRun::usage = "Helper for display";
pruneRunByDisks::usage = "Helper for display";
pruneRunToTopChain::usage = "Helper for display";
(*graphToContactLines::usage = "Helper for graphics";*)
diskNumbersInRun::usage = "Numbered (bare by default) disks in run";
bareNumberQ::usage = "False for left[] or right [] disks";
bareNumber::usage = "Take off left[] or right []";
getDiskFromRun::usage = "";
diskAndVisibleCopies::usage = "";
runDisks::usage = "Association of Disk[]s, by default ";
runDisksRadius::usage = "Association of radii of each disk";
runDisksHeight::usage = "Association of height of each disk";
runDisksDivergenceHeight::usage = "Association of xz of each disk";
diskR::usage = "better to provide an api for the functions that call this";
diskZ::usage = "";
diskXZ::usage = "";
postRunExtractNonOverlappingChains::usage = "";
deCycleChain::usage = "Used for parastichy counts and exposed as graphic helper"

(* debugging *)
influenceIntervals;
lowestByIntervals;
diskDiskIntersection;
diskdiskIntersectionQ;
globalRun::usage = "global _in this package_; made visible for debugging"; 


Begin["Private`"]


(* ::Input:: *)
(**)


(* ::Section:: *)
(*runmaking*)


(*runFromRun[run_] :=  Module[{g,d,chain},
	chain = Last[run["RunChains"]];
	g= chain["Chain"];
	d = Map[ <|#->run["DiskData"][#]|>&,Select[VertexList[g],bareNumberQ]];
	<|"Arena"->run["Arena"],"ContactGraph"->g,"DiskData"->d|>
];
*)
runFromLattice[lattice_] :=  Module[{g,d},
	If[!(lattice["d"]==0 && lattice["h"] == 1),
	Return[graphFromTCLattice[lattice]]];
	g= Graph[ {left[1]\[UndirectedEdge]1,1\[UndirectedEdge]right[1]}];
	d = <|1-><|"DiskNumber"->1,"Disk"->Disk[{0,0},0.5]|>|>;
	<|"ContactGraph"->g,"DiskData"->d|>
];

graphFromTCLattice[lattice_] := Module[{m,n,r,disks,diskx,diskxy,disknxy,disksOn,disksOf,diskData,g},
	{m,n}=First@Keys[lattice["parastichyNumbers"]];
	r=Sqrt[N@lattice["parastichyNumbers"][{m,n}]] ;
	disks= Map[Disk[#,r/2]&,lattice["namedLatticePoints"]];
	diskx[Disk[{x_,_},_]]:=x;
	diskxy[Disk[xy_,_]]:=xy;
	disknxy[node_] := diskxy[disks[node]];
	disknxy[left[node_]] := diskxy[disks[node]]+{-1,0};
	disknxy[right[node_]] := diskxy[disks[node]]+{1,0};
	disksOn= Flatten[Map[{#\[DirectedEdge]#-m,#\[DirectedEdge]#-n}&,Keys[lattice["namedLatticePoints"]]]];
	disksOf[a_\[DirectedEdge]b_] := Module[ {d1,d2},
		If[b<0,Return[Nothing[]]];
		d1=disks[a];d2=disks[b];
		If[Abs[diskx[d2]-diskx[d1]] <1/2, Return[a\[UndirectedEdge]b]];
		If[diskx[d1]> 0,Return[a\[UndirectedEdge]right[b]], Return[a\[UndirectedEdge]left[b]]]
	];
	g=Graph[Map[disksOf,disksOn]];
	g= Graph[g,VertexCoordinates->Map[disknxy[#]&,VertexList[g]],VertexLabels->"Name"];
	diskData = Association@KeyValueMap[#1-><|"DiskNumber"->#1,"Disk"->#2|>&,disks];
	<|"ContactGraph"->g,"DiskData"->diskData|>
];

restartRunFromRun[run_] := Module[{res,lastDisk},
	res=pruneRunToTopChain[run];
	res = KeyDrop[res,{"TopChain","RunChains","SpecifiedDiskMax"}];
	res
]


(* ::Input::Initialization:: *)
runDisks[run_,bareOnly_:False] := Association@Map[#->getDiskFromRun[run,#]&,diskNumbersInRun[run,bareOnly]];
runDisksRadius[run_,bareOnly_:False] := Map[diskR,runDisks[run,bareOnly]];
runDisksHeight[run_,bareOnly_:False] := Map[diskZ,runDisks[run,bareOnly]];
runDisksDivergenceHeight[run_,bareOnly_:False] :=  Map[diskXZ,runDisks[run,bareOnly]];


(* ::Input::Initialization:: *)



linearInterpolatorByEndPoints[{hStart_,hEnd_},{rStart_,rEnd_}] := Block[{rSlope},
rSlope=(rEnd-rStart)/(hEnd-hStart);
Function[{h}, rStart+ Piecewise[ {
{0,h< hStart}
,{(h-hStart) * rSlope , h< hEnd}
, { (hEnd-hStart) * rSlope ,True}
}]]
];



(* ::Section:: *)
(*Run  pruning*)


(* for display convenience *)





pruneRun[run_,Zrange_] := Module[{res,nodes},
	res=run;
	nodes=leftAndRightNumbers@nodesToPruneTo[res,Zrange];
	pruneRunByNodes[run,nodes]
];
nodesToPruneTo[run_,Zrange_] := Module[{g,nodes},
	g=run["ContactGraph"];
	nodes=Select[VertexList[g],IntervalMemberQ[Interval[Zrange],AnnotationValue[{g,#},VertexCoordinates][[2]]]&];
	nodes
];
pruneRunByNodes[run_,nodes_] := Module[{res},
	res=run;res["ContactGraph"] = Subgraph[res["ContactGraph"] ,nodes];
	res["RunChains"]= Select[res["RunChains"],useChainQ[#Chain,nodes]&];
	res["DiskData"]= KeyTake[res["DiskData"],nodes];
	If[KeyMemberQ[res,"NodeStatistics"],
		res["NodeStatistics"] = KeySelect[run["NodeStatistics"],MemberQ[nodes,#]&]];
	res["Arena"]["CylinderLU"] = MinMax[diskZ[#Disk]&/@ res["DiskData"]];
res
];
useChainQ[chain_,nodes_] := IntersectingQ[VertexList[chain],nodes];


pruneRunByDisks[run_,diskNumbers_] := Module[{nodes,res},
nodes=Flatten[leftAndRightNumbers/@Range[diskNumbers]]; 
pruneRunByNodes[run,nodes]
];

pruneRunToTopChain[run_] := Module[{res,chain,nodes},
chain=Last[run["RunChains"]]["Chain"];nodes=Flatten[leftAndRightNumbers[VertexList[chain]]];pruneRunByNodes[run,nodes]
];



(* ::Section:: *)
(*Chain code*)


(* ::Subsection:: *)
(*Utilities *)


(* ::Input::Initialization:: *)
diskZ[Disk[{_,z_},_]] := z;
diskX[Disk[{x_,_},_]] := x;
diskXZ[Disk[{x_,z_},_]] := {x,z};
diskRightX[Disk[{x_,_},r_]] := x+r;
diskLeftX[Disk[{x_,_},r_]] := x-r;
diskTopZ[Disk[{_,z_},r_]] := z+r;
diskBottomZ[Disk[{_,z_},r_]] := z-r;
diskR[Disk[{_,_},r_]] := r;


moveNumberRight[n_] := right[n];
moveNumberRight[left[n_]] := n;
moveNumberLeft[n_] := left[n];
moveNumberLeft[right[n_]] := n;
leftAndRightNumbers[n_List] := Join@@{n,moveNumberRight/@n,moveNumberLeft/@n};




bareNumber[n_] := n /; bareNumberQ[n];
bareNumber[n_] := First[n] /; rightNumberQ[n];
bareNumber[n_] := First[n] /; leftNumberQ[n];

bareNumberQ[n_] := IntegerQ[n];
leftNumberQ[n_] := StringEndsQ[ToString[Head[n]],"left"];
rightNumberQ[n_] := StringEndsQ[ToString[Head[n]],"right"];


moveDiskRight[Disk[{x_,z_},r_]] := Disk[{x+1,z},r];
moveDiskLeft[Disk[{x_,z_},r_]] := Disk[{x-1,z},r];
moveNumberedDiskRight[n_->d_] := moveNumberRight[n]->moveDiskRight[d];
moveNumberedDiskLeft[n_->d_] := moveNumberLeft[n]->moveDiskLeft[d];


diskNumbersInRun[run_,bareOnly_:True] := If[bareOnly,
	Keys[run["DiskData"]]
	,
	VertexList[run["ContactGraph"]]
	];

getDisk[n_] := getDiskFromRun[globalRun,n];
getDiskFromRun[run_,n_] := run["DiskData"][n]["Disk"] /; bareNumberQ[n];
getDiskFromRun[run_,n_] :=  moveDiskRight[getDiskFromRun[run,bareNumber[n]]] /; rightNumberQ[n];
getDiskFromRun[run_,n_] := moveDiskLeft[getDiskFromRun[run,bareNumber[n]]] /; leftNumberQ[n];

getDiskXZ[run_,n_] := diskXZ[getDiskFromRun[run,n]];
getDiskXZ[n_] :=getDiskXZ[globalRun,n] ;
getDiskZ[run_,n_] := diskZ[getDiskFromRun[run,n]];getDiskZ[n_] :=getDiskZ[globalRun,n] ;

diskIsLeft[Disk[{x_,z_},r_]] :=x<-1/2;
diskIsRight[Disk[{x_,z_},r_]] :=x>1/2;
diskIsNormalised[d_] := (! diskIsLeft[d]) && (! diskIsRight[d]);


disksMaximum[] := 
Max@Map[diskTopZ[getDisk[#]]&,Keys[globalRun["DiskData"]]];
diskHighestBottom[] := 
Max@Map[diskBottomZ[getDisk[#]]&,Keys[globalRun["DiskData"]]];

highestDiskZ[run_] := Max@Map[diskZ[getDiskFromRun[run,#]]&,Keys[run["DiskData"]]];

nextRadius[] := Module[{highestZ},
	highestZ=highestDiskZ[globalRun];
	globalRun["Arena"]["rFunction"][highestZ]
];
nextDiskNumber[] := Max[bareNumber/@ VertexList[globalRun["ContactGraph"]]]+1;


(* ::Section:: *)
(*run set up*)


executeRun[run_] := Module[{i,imax,res,diskTime,chain},
	Off[SSSTriangle::tri];
	globalRun=run;
		 
	globalRun["SpecifiedDiskMax"]= Max[bareNumber/@ VertexList[globalRun["ContactGraph"]]];
	globalRun["TopChain"] = topChainInRun[globalRun,globalRun["ContactGraph"]];

	globalRun["CurrentDiskChain"]=topChainToDiskChain[globalRun,globalRun["TopChain"]];
	globalRun["RunChains"]= Association[]; 
	
	If[!NumericQ[imax],imax=20000];
	
	monitorFunction := For[i=1,i<= imax,i++,
		If[runCompletesArena[],Break[]];
		diskTime= First@Timing[addNextDisk[]];
	];

	monitorString := Module[{z},
		z=diskHighestBottom[];
		StringTemplate["Next disk: `disk`;\nParastichy: `parastichy`;\nZ left: `zToMax`;\nr: `radius`;\nper-disk time: `timing`"][<|
			"disk"-> nextDiskNumber[]
			,"parastichy"-> reportLatestParastichy[]
			,"zToMax"->run["Arena"]["zMax"]-z
			,"radius"->diskR[Last[globalRun["DiskData"]]["Disk"]]
			,"timing"->diskTime
			|>]
	];

	If[False,monitorFunction,Monitor[monitorFunction,monitorString]];
	
	res= postRun[globalRun];
	res
];





runCompletesArena[] := Module[{cutoff},
	diskHighestBottom[]>globalRun["Arena"]["zMax"] || 
	Max[Keys[globalRun["DiskData"]]]> globalRun["Arena"]["diskMax"]
];




(* ::Section:: *)
(*addNextDisk*)


diskListRowFromDisk[node_,disk_] := Module[{},
node-><|"DiskNumber"->node,"Disk"->disk|>
];

addNextDisk[]:=Module[{nextR,row,n,timing,disk,tw,twc,xnextChain},
	nextR=nextRadius[];
	row=findNextDiskFromDiskChain[nextR];
	
	disk = row["NextDisk"];
	disk= jiggleDisk[globalRun,disk];
	n=nextDiskNumber[];
	updateContactGraph[n,row["NextDiskRestsOn"]];
	AppendTo[globalRun["DiskData"],diskListRowFromDisk[n,disk]];
	methods=globalRun["Arena"]["Methods"];
	
	nextdiskChain=modifyDiskChain[globalRun["CurrentDiskChain"],n,row]; 
	globalRun["CurrentDiskChain"]=nextdiskChain;
	
	logSupportChain[n];

	monitorZ= diskHighestBottom[];
	monitorRise=disksMaximum[]-monitorZ;

];

logSupportChain[n_] := Module[{chain,logEntry},
	chain=diskChainToGraphChain[globalRun["CurrentDiskChain"]];
	
	logEntry=<|
		"Chain"->chain,
		"Parastichy"->graphChainParastichy[chain]
		|>;

	AppendTo[globalRun["RunChains"],n->logEntry ];
	
	];
	
	diskChainParastichy


modifyDiskChain[diskChain_,n_,nextDiskRow_] := Module[{res,sortedSupportDiskNumbers,supportDisks},
	chainNumbers=Keys[diskChain];
	supportDisks= AssociationThread[ nextDiskRow["NextDiskRestsOn"]-> {
	nextDiskRow["Disk1"],nextDiskRow["Disk2"]}];
	supportDisks=SortBy[supportDisks,diskXZ[#][[1]]&];
	sortedSupportDiskNumbers=Keys[supportDisks];
	lrNumbers=Join[left/@chainNumbers,chainNumbers,right/@chainNumbers];
	{nl,nr}=sortedSupportDiskNumbers;
	
	toL= First@First@Position[lrNumbers,nl,1];
	fromR= First@First@Position[lrNumbers,nr,1];
	lrNumbers=Join[Take[lrNumbers,toL],{n},Drop[lrNumbers,fromR-1]];
	lookup=Append[diskChain,n-> nextDiskRow["NextDisk"]];
	lookup=KeyTake[lookup,lrNumbers];
	lookup
]


prettyGraph[run_,gp_] := Module[{g,setGraphXY,res,vc,es},
	g=gp;
	vc = Map[getDiskXZ[run,#]&,VertexList[g]];
	g = Graph[g,VertexCoordinates->vc];
	es = Map[#->edgeStyle[run,#]&,EdgeList[g]];
	g= Graph[g,EdgeStyle->es];
	g=Graph[g,VertexLabels->"Name"];
	g=Graph[g,VertexSize->Tiny];
	g
];

topChainToDiskChain[run_,chain_] := Module[{diskChain,lNode,rNode,path},
	diskChain= Map[#->getDiskFromRun[run,#]&,VertexList[chain]];	
	diskChain= KeySelect[diskChain,bareNumberQ];
	diskChain=SortBy[diskChain,#[[1,1]]&];
	lNode=First[Keys[diskChain]];
	rNode=Last[Keys[diskChain]];
	path=First@FindPath[chain,lNode,rNode];
	diskChain=KeyTake[diskChain,path];
	diskChain
];
diskChainToGraphChain[diskChain_] := Module[{lookup,g,nodes,edgeStyle,edgeColour,xz,edges,lastNode,leftDisk},
	nodes=Keys[diskChain];
	edges=UndirectedEdge@@@Partition[nodes,2,1];
	lastNode=Last[nodes];
	leftDisk=moveDiskLeft[diskChain[lastNode]];
	lookup=Append[diskChain,left[lastNode]-> leftDisk];
	nodes=Append[nodes,left[lastNode]];
	edges=Append[edges,left[lastNode]\[UndirectedEdge] First[nodes]];
	g=Graph[nodes,edges,VertexLabels->"Name",VertexSize->Tiny];

	xz=Map[diskXZ,lookup];
	edgeColour[a_ \[UndirectedEdge] b_] := If[xz[a][[2]]<xz[b][[2]],Red,Blue];
	edgeStyle =Map[#->edgeColour[#]&,EdgeList[g]];
	AnnotationValue[g,VertexCoordinates]=Normal[xz];
	AnnotationValue[g,EdgeStyle]=edgeStyle;
	g
];


graphChainParastichy[chain_]:=Module[{},
	res=Counts[Map[AnnotationValue[{chain,#},EdgeStyle]&,EdgeList[chain]]];
	(* in case left or right is missing *)
	res = Append[<|Red->0,Blue->0|>,res];
	res=KeySort[res];
	res

];



jiggleDisk[run_,disk_] := Module[{res,noise},
res=disk;
noise=run["Arena"]["Noise"];
If[MissingQ[noise],Return[res]];

res=jiggleDiskRadius[res,noise];
Return[res];
];

jiggleDiskRadius[Disk[xy_,r_],noise_] := Disk[xy,RandomReal[r*{1-noise,1+noise}]]




(* ::Subsection:: *)
(*findNextDiskFromSupportSet*)


findNextDiskFromChain[nextR_]:= Module[{chainNumbers,supportTable,res},
	chainNumbers=VertexList[globalRun["TopChain"]];
	supportTable=supportPairsFromSupportChainNumbers[chainNumbers,nextR];
	supportTable= SortBy[supportTable,diskZ[#NextDisk]&];
	If[Length[supportTable]==0,Print["No valid supports looking for ", nextDiskNumber[]];
Abort[]];
	res=First[supportTable];
	res
];

findNextDiskFromDiskChain[nextR_]:= Module[{supportTable,res},
	supportTable=supportPairsFromDiskChain[globalRun,nextR];

	supportTable= SortBy[supportTable,diskZ[#NextDisk]&];
	If[Length[supportTable]==0,Print["No valid supports looking for ", nextDiskNumber[]];Abort[]];
	res=First[supportTable];
	res
];

supportPairsFromDiskChain[run_,nextR_] :=Module[{supportDisks,pairs,supportTable,supportLine,supportLineFunction,newdisksToCheckIntersection,diskAboveSupportLineQ,extendedDisks,res,i,j},
	supportChainNumbers=Keys[run["CurrentDiskChain"]];
	supportDisks =addLeftRightSupporters[supportChainNumbers,nextR]; 
	supportDisks = Map[#->getDiskFromRun[run,#]&,supportDisks];
	pairs = Flatten[Table[{supportDisks[[i]],supportDisks[[j]]},{i,1,Length[supportDisks]},{j,i+1,Length[supportDisks]}],1];
	supportTable = Map[newMakeSupportTableRow[#]&,pairs];
	supportTable = Select[supportTable,Abs[#GapSeparation]< 2* nextR&];
	supportTable = Map[Append[#,"NextDisk"->newDiskOnThisPair[{#Disk1,#Disk2},nextR]]&,supportTable];
	supportTable = Select[supportTable,!MissingQ[#NextDisk]&];
	supportTable = Select[supportTable,isLeftRightSupported];
	supportTable = Select[supportTable,diskIsNormalised[#NextDisk]&];
	newdisksToCheckIntersection=Association[supportDisks];
(* the centre of the new disk must be above the line through the centres of the support *) 
	supportLine =  SortBy[Values@Map[diskXZ,newdisksToCheckIntersection],First];
	supportLineFunction = Interpolation[supportLine,InterpolationOrder->1];
	diskAboveSupportLineQ[disk_]:= diskZ[disk] > supportLineFunction[diskX[disk]];
	supportTable = Select[supportTable,diskAboveSupportLineQ[#NextDisk]&];
	supportTable= Map[newsupportTableFindIntersections[#,nextR,newdisksToCheckIntersection]&,supportTable];
	supportTable= Select[supportTable,Length[#Intersection]==0&];
	supportTable
];

isLeftRightSupported[row_] := Module[{x1,x2,xDisk,z1,z2,zDisk},
x1=diskX[row["Disk1"]];
x2 = diskX[row["Disk2"]];
xDisk = diskX[row["NextDisk"]];
z1=diskZ[row["Disk1"]];
z2 = diskZ[row["Disk2"]];
zDisk = diskZ[row["NextDisk"]];

IntervalMemberQ[Interval[{x1,x2}],xDisk] && (zDisk > z1 || zDisk > z2)
];


addLeftRightSupporters[sDisks_,nextR_] := Module[{supportDisks},
	supportDisks=sDisks;
	supportDisks= DeleteDuplicates[bareNumber/@ supportDisks];
	supportDisks =Flatten[Map[leftRightCouldSupport[#,nextR]&,supportDisks]];
	supportDisks
];

leftRightCouldSupport[node_,nextR_] := Module[{d,xy,res,nodeRadius},
	res= {node};
	d=getDisk[node];
	xy=diskXZ[d];
	nodeRadius=diskR[d];
	If[xy[[1]]-1 >=-0.5-  (nextR+nodeRadius),res=Append[res,left[node]]];
	If[xy[[1]]+1 <=   0.5+  (nextR+nodeRadius),res=Append[res,right[node]]];
	res 
];



(* ::Subsection:: *)
(*topChain*)


(* ::Input::Initialization:: *)
vectorXZ[run_,n1_\[DirectedEdge]n2_ | n1_\[UndirectedEdge] n2_ | {n1_,n2_}] := Module[{}, getDiskXZ[run,n2]-getDiskXZ[run,n1]];


neighbours[g_,n_] :=  Complement[VertexList[Graph[EdgeList[g,n\[UndirectedEdge]_]]],{n}]

vectorsFromNode[run_,g_,n_] := Module[{nbrs},
nbrs := Join[Map[ n\[DirectedEdge] # & , neighbours[g,n]],Map[ left[n]\[DirectedEdge] # & , neighbours[g,left[n]]]];
Association@Map[#->vectorXZ[run,#]&,nbrs]
];

clockwiseSortFunction[xy_] := Module[{x,y},
{x,y}=xy/Norm[xy];
-If[x>0,1+y,-(1+y)]
]; (* sorts in clockwise order *)

mostClockwiseEdgeFromNode[run_,g_,n_] := Module[{vectors},
vectors = vectorsFromNode[run,g,n];
If[Length[vectors]==0,Return[Missing[]]];
vectors = SortBy[vectors,clockwiseSortFunction];
First[Keys[vectors]]
];
(*
topChain[graph_] := topChainNumbersInRun[globalRun,graph];
*)
topChainInRun[run_,graph_] := 
Subgraph[graph,topChainNumbersInRun[run,graph] ];


(* ::Input::Initialization:: *)
topChainNumbersInRun[run_,graph_] := Module[{g,chain,neighbours,highestNode,lastNode,nextNode,from,nextEdge,chainMax=500},

g=graph;
highestNode=Last[SortBy[Select[VertexList[g],bareNumberQ],getDiskZ[run,#]&]];


nbrs[n_] := Join[Map[ n\[DirectedEdge] # & , neighbours[g,n]],Map[ left[n]\[DirectedEdge] # & , neighbours[g,left[n]]]];
chain= {highestNode};
debug=highestNode==59;debug=False;
For[i=1,i<chainMax,i++,

lastNode=Last[chain];
nextEdge =mostClockwiseEdgeFromNode[run,g,lastNode];

If[debug,Print[chain,nextEdge,prettyGraph[globalRun,g]]];If[MissingQ[nextEdge], (* may only happen for 1 left[1] megabodge *)
Break[]];
{from,nextNode}=List@@nextEdge;
If[MatchQ[from,left[_]],chain=Append[chain,from]];
chain=Append[chain,nextNode];
g=EdgeDelete[g,from\[UndirectedEdge]nextNode]; (* could move up into backwards branch *);
If[MatchQ[Last[chain],right[_]],
chain=Append[chain,moveNumberLeft[Last[chain]]]];

If[Last[chain]==First[chain],Break[];]
];If[i==chainMax,Print["i10"];Abort[]];


chain
];




(* ::Input:: *)
(**)


(* ::Subsection:: *)
(*supportPairsFromTopChain*)


reportLatestParastichy[] := Module[{chain},
chain = Values[globalRun["RunChains"]];
If[Length[chain]==0,Return[""]];
chain = Last[chain]["Parastichy"];
chain 
];




(* ::Input::Initialization:: *)
supportPairsFromSupportChainNumbers[supportChainNumbers_,nextR_] :=Module[{supportDisks,pairs,supportTable,supportLine,supportLineFunction,newdisksToCheckIntersection,diskAboveSupportLineQ,extendedDisks,res,i,j},

supportDisks =addLeftRightSupporters[DeleteDuplicates[supportChainNumbers],nextR]; 

supportDisks = Map[#->getDiskFromRun[globalRun,#]&,supportDisks];
pairs = Flatten[Table[{supportDisks[[i]],supportDisks[[j]]},{i,1,Length[supportDisks]},{j,i+1,Length[supportDisks]}],1];

supportTable = Map[newMakeSupportTableRow[#]&,pairs];
supportTable = Select[supportTable,Abs[#GapSeparation]< 2* nextR&];
supportTable = Map[Append[#,"NextDisk"->newDiskOnThisPair[{#Disk1,#Disk2},nextR]]&,supportTable];
supportTable = Select[supportTable,!MissingQ[#NextDisk]&];
supportTable = Select[supportTable,isLeftRightSupported];
supportTable = Select[supportTable,diskIsNormalised[#NextDisk]&];


newdisksToCheckIntersection=Association[supportDisks];
(* the centre of the new disk must be above the line through the centres of the support *) 

supportLine =  SortBy[Values@Map[diskXZ,newdisksToCheckIntersection],First];
supportLineFunction = Interpolation[supportLine,InterpolationOrder->1];
diskAboveSupportLineQ[disk_]:= diskZ[disk] > supportLineFunction[diskX[disk]];

supportTable = Select[supportTable,diskAboveSupportLineQ[#NextDisk]&];


supportTable= Map[newsupportTableFindIntersections[#,nextR,newdisksToCheckIntersection]&,supportTable];
supportTable= Select[supportTable,Length[#Intersection]==0&];

supportTable
];

isLeftRightSupported[row_] := Module[{x1,x2,xDisk,z1,z2,zDisk},
x1=diskX[row["Disk1"]];
x2 = diskX[row["Disk2"]];
xDisk = diskX[row["NextDisk"]];
z1=diskZ[row["Disk1"]];
z2 = diskZ[row["Disk2"]];
zDisk = diskZ[row["NextDisk"]];

IntervalMemberQ[Interval[{x1,x2}],xDisk] && (zDisk > z1 || zDisk > z2)
];


newMakeSupportTableRow[{node1_->d1_,node2_->d2_}] := Module[{res,hDiff,eSep},
hDiff= diskX[d2]-diskX[d1];

If[hDiff> 0,
eSep= diskLeftX[d2]-diskRightX[d1] 
,
eSep = diskLeftX[d1]-diskRightX[d2] 
];
eSep = If[eSep>0,eSep,0];
<|"Disk1"->d1,"Disk2"->d2,"NextDiskRestsOn"->{node1,node2}, "GapSeparation"->eSep|>
]

newsupportTableFindIntersections[row_,nextR_,disks_] := Module[{res,extendedDisks,locationIntersectsQ,intersections},
intersections= Map[diskdiskIntersectionQ[row["NextDisk"],#]&,KeyDrop[disks,row["NextDiskRestsOn"]]];
intersections= Keys@Select[ intersections,TrueQ];
Append[row,"Intersection"->intersections]
];


newDiskOnThisPair[{d1_,d2_},r_] := Module[{res},
res = diskdiskUpperTouchingPoint[{d1,d2},r];
If[MissingQ[res],Return[res]]; (* no overlap *) 
Disk[res,r]
];

diskdiskIntersectionQ[Disk[xy1_,r1_],Disk[xy2_,r2_]]:= Norm[xy1-xy2,2] < (r1+ r2);


disksInSupportChain[chain_,nextR_] :=
addLeftRightSupporters[DeleteDuplicates[chain],nextR]; 



(* ::Input::Initialization:: *)
isLeftRightSupported[row_] := Module[{x1,x2,xDisk,z1,z2,zDisk},
x1=diskX[row["Disk1"]];
x2 = diskX[row["Disk2"]];
xDisk = diskX[row["NextDisk"]];
z1=diskZ[row["Disk1"]];
z2 = diskZ[row["Disk2"]];
zDisk = diskZ[row["NextDisk"]];

IntervalMemberQ[Interval[{x1,x2}],xDisk] && (zDisk > z1 || zDisk > z2)
];


newMakeSupportTableRow[{node1_->d1_,node2_->d2_}] := Module[{res,hDiff,eSep},
hDiff= diskX[d2]-diskX[d1];

If[hDiff> 0,
eSep= diskLeftX[d2]-diskRightX[d1] 
,
eSep = diskLeftX[d1]-diskRightX[d2] 
];
eSep = If[eSep>0,eSep,0];
<|"Disk1"->d1,"Disk2"->d2,"NextDiskRestsOn"->{node1,node2}, "GapSeparation"->eSep|>
]

newsupportTableFindIntersections[row_,nextR_,disks_] := Module[{res,extendedDisks,locationIntersectsQ,intersections},
intersections= Map[diskdiskIntersectionQ[row["NextDisk"],#]&,KeyDrop[disks,row["NextDiskRestsOn"]]];
intersections= Keys@Select[ intersections,TrueQ];
Append[row,"Intersection"->intersections]
];


newDiskOnThisPair[{d1_,d2_},r_] := Module[{res},
res = diskdiskUpperTouchingPoint[{d1,d2},r];
If[MissingQ[res],Return[res]]; (* no overlap *) 
Disk[res,r]
];

diskdiskIntersectionQ[Disk[xy1_,r1_,z___],Disk[xy2_,r2_,z___]]:= Norm[xy1-xy2,2] < (r1+ r2);




(* ::Subsection:: *)
(*updateContactGraph*)


updateContactGraph[n_,{n1_,n2_}] := Module[{g,d,d1,d2,d1Left},
g = globalRun["ContactGraph"];
g= updateGraph[g,n,{n1,n2}];
globalRun["ContactGraph"]=g;
];

SetAttributes[updateGraph,HoldFirst];
updateGraph[g_,n_,{n1_,n2_}] := Module[{},
g=  VertexAdd[g,n];
If[!bareNumberQ[n1],
 g= VertexAdd[g,n1]];
If[!bareNumberQ[n2],
 g= VertexAdd[g,n2]];
g= EdgeAdd[g,
{n \[UndirectedEdge] n1,n \[UndirectedEdge] n2}];
g
];



(* ::Section:: *)
(*Interval calculations*)


chainToIntervals[diskChain_,r_]  := Module[{inhibitionData},
	inhibitionData=KeyValueMap[
		<|
		"Interval"->Interval[{diskLeftX[#2]-r,diskRightX[#2]+r}],
		"Form"-><|#1->Disk[diskXZ[#2],r+diskR[#2]]|>
		|>&,diskChain];
	inhibitionData
];


lowestByIntervals[diskChain_,r_,n_:All]:= Module[{intervalList,lowest,lowestLeftPosition,lowestLeftDiskNumber,lowestRightPosition,lowestRightDiskNumber},
	intervalList=influenceIntervals[diskChain,r,n];
	(*Print@frontDisplayer[diskChain,intervalList];
*)
	lowest=First@SortBy[intervalList,#RightZ&];
	lowestLeftPosition=First@First@Position[intervalList,lowest];
	lowestLeftDiskNumber=First@Keys[intervalList[[lowestLeftPosition]]["Form"]];
	If[lowestLeftPosition==Length[intervalList],
		lowestRightPosition=1
		,
		lowestRightPosition=lowestLeftPosition+1
		];
	lowestRightDiskNumber=First@Keys[intervalList[[lowestRightPosition]]["Form"]];
	<|"SupportPair"->{lowestLeftDiskNumber,lowestRightDiskNumber}|>
		
];



influenceIntervals[diskChain_,r_,n_:All]:= Module[{inhibitionData,intervalList,res,minZ},
	inhibitionData=chainToIntervals[diskChain,r];
	
	forms=Map[First[Values[#Form]]&,inhibitionData];
	minZ=Min@Map[diskBottomZ,forms];
	intervalList= {
	<|"Interval"->Interval[{-1/2,1/2}],
	"Form"-><|-\[Infinity]->Line[{{-1/2,minZ},{1/2,minZ}}]|>|>
	};
	
	
	intervalList=Fold[addResolvedInfluences,intervalList,
		Take[inhibitionData,n]];
	intervalList=Map[KeySort,intervalList];
	
	intervalList=mergeIntervalList[intervalList];
	intervalList=Map[Append[#,
	<|"RightZ"->formAtIntervalEnd[#]	|>]&,intervalList];
	intervalList
];
formAtIntervalEnd[interval_] := Module[{x,z,disk},
	disk=First@Values@interval["Form"];
	x= Max[interval["Interval"]];
	If[x==1/2 && !(diskX[disk]+diskR[disk]==1/2),
		Return[\[Infinity]]];
	z= diskZAtX[disk,x];
	z
];


adjacentFormQ[<|___,"Form"->form1_,___|>,<|___,"Form"->form2_,___|>] := Module[{},First[Keys[form1[]]]==First[Keys[form2[]]]];
mergeAdjacentForms[<|___,"Form"->form1_,"Interval"->interval1_|>,<|___,"Form"->form2_,"Interval"->interval2_|>] := 
	<|"Form"->form1,"Interval"->Interval[{Min[interval1],Max[interval2]}]|>;
mergeIntervalList[intervalList_] := Module[{res},
	res=intervalList//. ({b___,f1_,f2_,c___}:>{b,mergeAdjacentForms[f1,f2],c} /; adjacentFormQ[f1,f2]);
	res
];


addResolvedInfluences[intervalList_,influence_] := Module[{res},

	multiFormIntervals=addInfluencesToIntervals[intervalList,influence];
	res= Flatten@Map[resolveFormsToForm,multiFormIntervals];
	
	res 
	];

frontDisplayer[diskChain_,intervalList_] := (ix=0;Graphics[
	{
	{FaceForm[LightGray],Values[diskChain],
	KeyValueMap[Text[#1,#2]&,diskXZ/@diskChain],
	InfiniteLine[{-1/2,0},{0,1}],InfiniteLine[{1/2,0},{0,1}]},
	Map[formDisplayer,intervalList]},
	Axes->True,PlotRange->All]);
	
formDisplayer[interval_] := Module[{},
		{
		ColorData[17][ix++],PointSize[Large],
		Point[{Max[interval["Interval"]],interval["RightZ"]}],
	InfiniteLine[{Max[interval["Interval"]],0},{0,1}],
	Select[Values[interval["Form"]],Head[#]=!= Line &]/. Disk->Circle
	}
];

emptyIntervalQ[Interval[interval_]]:= False;
emptyIntervalQ[Interval[interval___]]:= True;
	
addInfluencesToIntervals[intervals_,influence_]:= Module[{res},
	res=Flatten@Map[addInfluenceToInterval[#,influence]&,intervals];
	res
];


addInfluenceToInterval[interval_,influence_]:= Module[
	{intervalLeft,intervalRight,influenceLeft,influenceRight,leftInterval,middleInterval,rightInterval,intervals},
	debug=First@Keys@influence["Form"]==56;debug=False;
	If[debug,Print["\n\nInterval: ",interval,"\nInfluence: ",influence]];
	If[emptyIntervalQ[IntervalIntersection[interval["Interval"],influence["Interval"]]],If[debug,Print["Empty"]];
		Return[interval]];
	{intervalLeft,intervalRight}=MinMax[interval["Interval"]];
	{influenceLeft,influenceRight}=MinMax[influence["Interval"]];
	leftInterval =If[influenceLeft<intervalLeft,
		Missing[]
		,
		<|"Form"-><|interval["Form"]|>,"Interval"->Interval[{intervalLeft,influenceLeft}]|>
		];
	middleInterval=
		<|"Form"->Join[interval["Form"],influence["Form"]],
		"Interval"->Interval[{Max[intervalLeft,influenceLeft],Min[intervalRight,influenceRight]}]|>;
	rightInterval= If[influenceRight>intervalRight,
		Missing[]
		,
		<|"Form"-><|interval["Form"]|>,"Interval"->Interval[{influenceRight,intervalRight}]|>
		];
	If[debug,Print["InfluenceIntervals:",{leftInterval,middleInterval,rightInterval}]];
	intervals=DeleteMissing[{leftInterval,middleInterval,rightInterval}];

	intervals
];


resolveFormsToForm[interval_] := Module[{forms,res},
	
	If[Length[interval["Form"]]==1,
			Return[interval]
	];
	If[Length[interval["Form"]]>2, Print["Multiform abort"];Abort[]];

	res=resolveTwoForms[interval];
	res
	]


resolveTwoForms[interval_] := Module[{res,disk1,disk2,nonIntersectingIntervals,intersectionX,intersectionXZ},
	(* intervals on which two disks don't intersect are separated by where they do intersect *)
	{disk1,disk2}=Values[interval["Form"]];
	debug=Keys@interval["Form"]=={46,56};debug=False;
	If[debug,Print["Resolving ",interval]];
	If[Head[disk1]===Line,
		res=interval;
		res["Form"]=KeyDrop[interval["Form"],{Key[-\[Infinity]]}];
		Return[res];
		];
	intersectionXZ=diskDiskIntersection[disk1,disk2,interval];
	intersectionX=First/@intersectionXZ;
	intersectionX= Sort@Join[{
		Min[interval["Interval"]]},
		intersectionX,
		{Max[interval["Interval"]]}
	];
	If[debug,Print["iXZ",intersectionXZ,intersectionX]];
	nonIntersectingIntervals=Interval/@Partition[intersectionX,2,1];
	nonIntersectingIntervals=Map[<|"Interval"->#,"Form"->interval["Form"]|>&,nonIntersectingIntervals];
	nonIntersectingIntervals=Map[resolveNonIntersectingForms,nonIntersectingIntervals];
	If[debug,Print[" to ",nonIntersectingIntervals]];
	Return[nonIntersectingIntervals];
	
	]


resolveNonIntersectingForms[interval_] := Module[{xl,xr,xmid,res,activeDiskNumber,numberedFormsZnumberedFormsZ},
	{xl,xr}=MinMax[interval["Interval"]];
	xmid=Mean[{xl,xr}];
	numberedFormsZ=Map[diskZAtX[#,xmid]&,interval["Form"]];
	numberedFormsZ=Reverse@Sort[numberedFormsZ];
	activeDiskNumber=First[Keys[numberedFormsZ]];
	res=interval;
	res["Form"]=KeyTake[res["Form"],{activeDiskNumber}];
	res
	];


intervalDisplay[intervals_] := Module[{},
	Print[intervals]
];

diskZAtX[Disk[{x1_,z1_},r1_],x_]:= Module[{det},
	det = Chop[r1^2-(x-x1)^2];
	res = z1+Sqrt[det];
	res
	];
diskZAtX[Line[{{a_,b_},{c_,d_}}],x_]:=b;
	
diskDiskIntersection[Disk[{x1_,y1_},r1_],Disk[{x2_,y2_},r2_],interval_]:= Module[{},
	If[Norm[{x1,y1}-{x2,y2},2] > (r1+ r2),
		Return[{}]
	];
	ix=circleCircleIntersection[Circle[{x1,y1},r1],Circle[{x2,y2},r2]];
	ix=Select[ix,IntervalMemberQ[(interval["Interval"]),First[#]]&];
	ix  
	
	]


circleCircleIntersection[Circle[{x1_,y1_},r1_],Circle[{x2_,y2_},r2_]] := Module[{sols,x,y},
	sols =
		{{x->(r2^2 (x1-x2)+r1^2 (-x1+x2)+(x1+x2) ((x1-x2)^2+(y1-y2)^2)-\[Sqrt](-(((r1-r2)^2-(x1-x2)^2-(y1-y2)^2) ((r1+r2)^2-(x1-x2)^2-(y1-y2)^2) (y1-y2)^2)))/(2 ((x1-x2)^2+(y1-y2)^2)),y->(-r1^2 (y1-y2)^2+x1 \[Sqrt](-((r1^4+(-r2^2+(x1-x2)^2+(y1-y2)^2)^2-2 r1^2 (r2^2+(x1-x2)^2+(y1-y2)^2)) (y1-y2)^2))-x2 \[Sqrt](-((r1^4+(-r2^2+(x1-x2)^2+(y1-y2)^2)^2-2 r1^2 (r2^2+(x1-x2)^2+(y1-y2)^2)) (y1-y2)^2))+(y1-y2) (r2^2 (y1-y2)+((x1-x2)^2+(y1-y2)^2) (y1+y2)))/(2 ((x1-x2)^2+(y1-y2)^2) (y1-y2))},{x->(r2^2 (x1-x2)+r1^2 (-x1+x2)+(x1+x2) ((x1-x2)^2+(y1-y2)^2)+\[Sqrt](-(((r1-r2)^2-(x1-x2)^2-(y1-y2)^2) ((r1+r2)^2-(x1-x2)^2-(y1-y2)^2) (y1-y2)^2)))/(2 ((x1-x2)^2+(y1-y2)^2)),y->(-r1^2 (y1-y2)^2-x1 \[Sqrt](-((r1^4+(-r2^2+(x1-x2)^2+(y1-y2)^2)^2-2 r1^2 (r2^2+(x1-x2)^2+(y1-y2)^2)) (y1-y2)^2))+x2 \[Sqrt](-((r1^4+(-r2^2+(x1-x2)^2+(y1-y2)^2)^2-2 r1^2 (r2^2+(x1-x2)^2+(y1-y2)^2)) (y1-y2)^2))+(y1-y2) (r2^2 (y1-y2)+((x1-x2)^2+(y1-y2)^2) (y1+y2)))/(2 ((x1-x2)^2+(y1-y2)^2) (y1-y2))}};
	SortBy[{x,y} /. sols ,Last]
];





(*iComplement=ResourceFunction["IntervalComplement"];

addInfluenceToIntervalType[interval_,disk_,{Line,Disk}] :=Module[{},
	unchangedIntervals=iComplement[
		interval["Interval"]
		,Interval[{diskLeftX[disk],diskRightX[disk]}]
	];
	unchangedIntervalSets= Interval/@ (List@@unchangedIntervals);
	unchangedIntervalSets=Map[<|"LeftZ"->interval["LeftZ"],
	"RightZ"->interval["RightZ"],"Interval"->#,"Form"->Line[{{Min[#],interval["LeftZ"]},{Max[#],interval["RightZ"]}}]|>&,unchangedIntervalSets];
	changedInterval=IntervalIntersection[
		interval["Interval"]
		,Interval[{diskLeftX[disk],diskRightX[disk]}]];
	If[diskLeftX[disk]<Min[changedInterval],
		diskIntervalLeftZ=diskEvaluate[disk,Min[changedInterval]];
		angleLeft=diskAngle[disk,{Min[changedInterval],diskIntervalLeftZ}]
	,
		diskIntervalLeftZ=diskZ[disk];
		angleLeft=\[Pi]
	];
	If[diskRightX[disk]>Max[changedInterval],
		diskIntervalRightZ=diskEvaluate[disk,Max[changedInterval]];
		angleRight=diskAngle[disk,{Max[changedInterval],diskIntervalRightZ}]
	,
		diskIntervalRightZ=diskZ[disk];
		angleRight=0
	];
	(* angleLeft etc only for display *)
	
	changedIntervalData={<|
	"LeftZ"->diskIntervalLeftZ,
	"RightZ"->diskIntervalRightZ,
	"Interval"->changedInterval,
	"Form"->Disk[diskXZ[disk],diskR[disk],{0,angleLeft}]
	|>};
	res=Join[unchangedIntervalSets,changedIntervalData];
	
		res
	];
	
showIntervals[intervals_] := Module[{zMin},
	zMin=Min@Map[KeyTake[#,{"LeftZ","RightZ"}]&,intervals];
	{LightGray,Line[{{-1/2,zMin},{1/2,zMin}}]
	,Map[#Form&,intervals]
	}
];*)

	


(*addInfluenceToIntervalType[interval_,disk_,{Disk,Disk}] :=Module[{unchangedIntervals},
	unchangedIntervals=iComplement[
		interval["Interval"]
		,Interval[{diskLeftX[disk],diskRightX[disk]}]
	];
	Print[unchangedIntervals];
	unchangedIntervalSets= Interval/@ (List@@unchangedIntervals);
	unchangedIntervalSets=Map[<|"LeftZ"->interval["LeftZ"],
	"RightZ"->interval["RightZ"],"Interval"->#,"Form"->Line[{{Min[#],interval["LeftZ"]},{Max[#],interval["RightZ"]}}]|>&,unchangedIntervalSets];
	
	
	changedInterval=IntervalIntersection[
		interval["Interval"]
		,Interval[{diskLeftX[disk],diskRightX[disk]}]];
	returnIntervalList ={
		interval;
	}	
		
	diskIntervalLeftZ=diskEvaluate[disk,Min[changedInterval]];
	diskIntervalRightZ=diskEvaluate[disk,Max[changedInterval]];
	
	
	If[diskIntervalLeftZ < interval["LeftZ"] && diskIntervalRightZ < interval["RightZ"],
		Print["New disk has no effect"]
		];

	If[diskLeftX[disk]<Min[changedInterval],
		diskIntervalLeftZ=diskEvaluate[disk,Min[changedInterval]];
		angleLeft=diskAngle[disk,{Min[changedInterval],diskIntervalLeftZ}]
	,
		diskIntervalLeftZ=diskZ[disk];
		angleLeft=\[Pi]
	];
	If[diskRightX[disk]>Max[changedInterval],
		diskIntervalRightZ=diskEvaluate[disk,Max[changedInterval]];
		angleRight=diskAngle[disk,{Max[changedInterval],diskIntervalRightZ}]
	,
		diskIntervalRightZ=diskZ[disk];
		angleRight=0
	];
	(* angleLeft etc only for display *)
	
	changedIntervalData= {<|
		"LeftZ"->diskIntervalLeftZ,
		"RightZ"->diskIntervalRightZ,
		"Interval"->changedInterval,
		"Form"->Disk[diskXZ[disk],diskR[disk],{0,angleLeft}]
	|>};
	Join[unchangedIntervalSets,changedIntervalData]
	];
		*)


(*addInfluenceToIntervalType[interval_,influence_,type_] := Module[{},
	Print[" type ", type, "unimplemented"];
	Print[interval,influence];
	Abort[]
	];
	
diskEvaluate[Disk[{x_,z_},r_,___],atX_] := Module[{},
	xDistance= Abs[x-atX];
	z+Sqrt[r^2-xDistance^2]
];

diskAngle[Disk[{x_,z_},r_,___],{atX_,atZ_}]  :=
 Module[{},
	xDistance= atX-x;
	zDistance=atZ-z;
	ArcTan[xDistance,zDistance]
]

	*)






(* ::Section:: *)
(*Disk calculations*)


(* ::Input::Initialization:: *)
(*  *)
diskdiskUpperTouchingPoint[pairDisks_,r_] := Module[
{lrPoints},
lrPoints= newdiskdiskTouchingPoint[pairDisks,r];

If[MissingQ[lrPoints],Return[lrPoints]];
Last[SortBy[lrPoints,N@Last[#]&]]
];



(* ::Input::Initialization:: *)



newdiskdiskTouchingPoint[diskPair_,r_] := Module[{c1,c2,r1,r2,interdisk,interdiskVectorNorm,interdiskNormal,
angle,vector,normal,interdiskVector},
{c1,c2}= Map[diskXZ,diskPair];
{r1,r2}= Map[diskR,diskPair];
interdiskVector = c2-c1;
interdisk = Norm[interdiskVector];
interdiskVectorNorm = interdiskVector/Norm[interdiskVector];
interdiskNormal =  {-interdiskVectorNorm[[2]],interdiskVectorNorm[[1]]};


angle = sssTriangleInteriorAngle[r2+r,r+r1,interdisk];
If[MissingQ[angle],Return[angle]];
vector=(r1+r)*Cos[angle]*interdiskVectorNorm;
normal = (r1+r)*Sin[angle]* interdiskNormal;
{c1+vector+normal,c1+vector-normal}
]

sssTriangleInteriorAngle[a_,b_,c_] := Module[{tri,angle},
	tri =SSSTriangle[a,b,c];
	If[Head[tri]==SSSTriangle,
	Return[Missing["Not a triangle"]]];
	angle = TriangleMeasurement[tri,{"InteriorAngle",1}];
	angle
];



(* ::Section:: *)
(*Post run prettification*)


postRun[run_] := Module[{res,chains,nrchains},
res = pretty[run];
edgeCheck[res];
(*res= postRunNodeStatistics[res];
*)
res
];

pretty[run_] := Module[{g,setGraphXY,res,vc,es},
	g=prettyGraph[run,run["ContactGraph"]];
	res=run;
	res["ContactGraph"]=g;
	res
];

edgeStyle[run_,upper_ \[UndirectedEdge] lower_] := Module[{vxy},
vxy=vectorXZ[run,upper\[UndirectedEdge] lower];
If[First[vxy]>= 0,Blue,Red]
]
edgeStyle[upper_ \[DirectedEdge] lower_] := edgeStyle[upper \[UndirectedEdge] lower];
edgeStyle[run_,upper_ \[DirectedEdge] lower_] := edgeStyle[run,upper \[UndirectedEdge] lower];


(* ::Input::Initialization:: *)



(* ::Input:: *)
(* edgeCheck[run_] := Module[{g,res,edgeStyler,res2,nodesToCheck,lopsidedAllowedTo},g=run["ContactGraph"];res=Map[lowerEdges[g,#]&,Complement[Select[VertexList[g],bareNumberQ],{1}]];If[Or@@ Map[Length[#]!=2&,res],"Print some nodes without two supports"];*)
(**)
(*edgeStyler[edgeList_] := Map[edgeStyle,edgeList];*)
(*lopsidedAllowedTo = run["SpecifiedDiskMax"];*)
(**)
(*nodesToCheck = Complement[Select[VertexList[g],bareNumberQ],Range[lopsidedAllowedTo]];res2= Map[lowerEdges[g,#]&,nodesToCheck];res2 =Association@Map[#->Sort@edgeStyler[#]&,res2];res2= Select[res2,# != {RGBColor[0, 0, 1],RGBColor[1, 0, 0]}&];*)
(**)
(*If[Length[res2]>0 &&  Keys[res2]!={{}},*)
(*Print["Print  lopsided nodes", res2];*)
(*];*)
(**)
(*];*)
(**)
(*edgeCheckNode[n_] := Module[{g,res,edgeStyler,res2,le},*)
(*If[bareNumber[n]===1,Return[True]];*)
(*g=globalRun["ContactGraph"];*)
(*le = lowerEdges[g,n];*)
(*If[Length[le]!=2,Print["Node ", n, " without two supports"]];*)
(*If[Sort[Map[edgeStyle,Keys[le]]]!= {RGBColor[0, 0, 1],RGBColor[1, 0, 0]},*)
(*Print["Node ", n, " is lopsided"]];*)
(**)
(**)
(*];*)
(**)
(*lowerEdges[g_,n_] := Module[{nxy,edges,nodes},*)
(*nxy=AnnotationValue[{g,n},VertexCoordinates];*)
(*edges= EdgeList[g,n \[UndirectedEdge]_];*)
(*nodes= Complement[Flatten[List@@@edges],{n}];*)
(*nodes= Association@Map[(n \[DirectedEdge] #)-> Last[AnnotationValue[{g,#},VertexCoordinates]]-Last[nxy]&,nodes];*)
(*nodes = Select[nodes,#<=0&];*)
(*Keys[nodes]*)
(*];*)


(* ::Input::Initialization:: *)
edgeDirectedRightwards[g_,v1_ \[UndirectedEdge] v2_] := Module[{v1xy,v2xy},
v1xy=AnnotationValue[{g,v1},VertexCoordinates];
v2xy=AnnotationValue[{g,v2},VertexCoordinates];
If[First[v1xy]<First[v2xy],
v1\[DirectedEdge] v2,v2\[DirectedEdge] v1]
];

directedGraphRightwards[g_] := Module[{v,vxy,edges,d},
v=VertexList[g];
vxy=Map[#->AnnotationValue[{g,#},VertexCoordinates]&,v];
edges =EdgeList[g];
edges= Map[edgeDirectedRightwards[g,#]&,edges];
d=Graph[v,edges,VertexCoordinates->vxy,VertexSize->Tiny,VertexLabels->"Name"
,EdgeShapeFunction->{{"HalfFilledArrow","ArrowSize"->.05}}
];
d
];


(* ::Input::Initialization:: *)
leftRightChainEdges[d_] := Module[{nonbares,nonbarePaths,nonbarePath,pathToDE,pathscore,chains},
nonbares=Select[VertexList[d],!bareNumberQ[#]&];
nonbarePaths[left[node_]] :=  FindPath[d,left[node],node,\[Infinity],All];
nonbarePaths[right[node_]] :=  FindPath[d,node,right[node],\[Infinity],All];
pathscore[path_] := Module[{z},
z= Map[Last@AnnotationValue[{d,#},VertexCoordinates]&,path];
z= z- Mean[z];
z= Map[Abs,z];
z= Total[z];
z
];
nonbarePath[node_] :=  First@SortBy[nonbarePaths[node],pathscore];

pathToDE[path_] := (monitorTopChainCount=path;Map[Apply[UndirectedEdge,#]&,Partition[path,2,1]]);
chains=Association@Map[(monitorFlatChain=#;bareNumber[#]->pathToDE[nonbarePath[#]])&,nonbares];
chains
];
(*
pickShortest[edgeList_] := First@SortBy[edgeList,#Xdifference&];

*)





(* ::Input::Initialization:: *)
countbyPath[chain_] := Module[{res},
res=Counts[classifyPath[chain]];
res= Append[<|RGBColor[0, 0, 1]->0,RGBColor[1, 0, 0]->0|>,res];
res=KeySort[res];
res
]; 
classifyPath[chain_] := Map[classifyPathEdge[chain,#]&,EdgeList[chain]];
classifyPathEdge[g_,v1_\[UndirectedEdge] v2_] := Module[{v1xy,v2xy},
v1xy=AnnotationValue[{g,v1},VertexCoordinates];
v2xy=AnnotationValue[{g,v2},VertexCoordinates];
{v1xy,v2xy}= SortBy[{v1xy,v2xy},First];
If[ Last[v1xy]>=Last[v2xy],Red,Blue]
];



(* Obsolete *)
(*
flattestChains[g_] := Module[{d,chainEdges,chains},
d=directedGraphRightwards[g];
chainEdges= leftRightChainEdges[d];
chains= Map[<|"Chain"->Subgraph[g,#]|>&,chainEdges];
chains = Map[Append[#,"Parastichy"->countbyPath[#Chain]]&,chains];
chains
]
*)


(* ::Section:: *)
(*Post run statistics*)


(* ::Input::Initialization:: *)
(*postRunNodeStatistics[run_] := Module[{res,angles},
res = run;
angles= Map[<|"Angle"->#|>&,computeNodeLowerAngles[run]];
res= Append[res,"NodeStatistics"-> angles];
res 
];*)


(* ::Section:: *)
(*Chains by flattest*)


(* ::Input::Initialization:: *)
(* Obsolete *)
(*
flattestChainData[graph_,chainsToDo_:100] := Module[{i,g,chainGraphSet,chains,chainSet},
chains= flattestChains[graph];

chainSet=chains /. DirectedEdge -> UndirectedEdge;
chainSet= Map[Append[#,"MeanZ"->chainMeanZ[#Chain]]&,chainSet];
chainSet= Map[Append[#,"MeanRadius"->chainMeanRadius[#Chain]]&,chainSet];
chainSet
];
*)


(*chainMeanZ[chain_] := Last@Mean@Map[AnnotationValue[{chain,#},VertexCoordinates]&,Drop[VertexList[chain],-1]];
*)
(*chainMeanRadius[chain_] := Module[{e,gxy,norm},
gxy[n_] :=AnnotationValue[{chain,n},VertexCoordinates];
norm[n1_,n2_] := Norm[gxy[n2]-gxy[n1]];
e=EdgeList[chain];
e=norm @@@ e;
Mean[e]/2
];*)



(* ::Section:: *)
(*Top down chains*)


(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
(*
chainMeanLatticeAngle[run_,chain_] := Module[{res,v},
v=Select[VertexList[chain],bareNumberQ];
res= KeyTake[run["NodeStatistics"],v];
res= Map[#Angle&,res];
Mean[res]
];
chainSDLatticeAngle[run_,chain_] := Module[{res,v},
v=Select[VertexList[chain],bareNumberQ];
res= KeyTake[run["NodeStatistics"],v];
res= Map[#Angle&,res];
If[Length[res]<2,Return[Missing["monoEdge"]]];
StandardDeviation[res]
];
*)
(*
chainMeanZ[chain_] := Last@Mean@Map[AnnotationValue[{chain,#},VertexCoordinates]&,Drop[VertexList[chain],-1]];

chainMeanRadius[chain_] := Module[{e,gxy},
gxy[n_] :=AnnotationValue[{chain,n},VertexCoordinates];
norm[n1_,n2_] := Norm[gxy[n2]-gxy[n1]];
e=EdgeList[chain];
e=norm @@@ e;
Mean[e]/2
];*)



(* ::Input::Initialization:: *)
(*stripChain[g_,chain_]:= Module[{lrchain},
lrchain =leftAndRightNumbers[DeleteDuplicates[bareNumber/@chain]];

lrchain= Intersection[lrchain,VertexList[g]];
VertexDelete[g,lrchain]
];
*)
chainEncircles[chain_] := Module[{lrNode},
lrNode=Select[VertexList[chain],!bareNumberQ[#]&];If[Length[lrNode]!=1,Return[False]];MemberQ[VertexList[chain],bareNumber[First[lrNode]]]
];

chainParastichy[chain_]:=Module[{aChain},
If[!chainEncircles[chain],Return[Missing["Not encircling"]]];

aChain = If[AcyclicGraphQ[chain],chain,deCycleChain[chain]];
If[MissingQ[aChain],Return[aChain]];
acyclicChainParastichy[aChain]

];

oldDeCycleChain[chain_] := Module[{reducedchain,v,first,last,edges,res},v=VertexList[chain];v=SortBy[v,AnnotationValue[{chain,#},VertexCoordinates]&];first=First[v];last=Last[v];
edges= Flatten[FindCycle[chain,Infinity,All],2];vInLoops=Complement[DeleteDuplicates[Flatten[List@@@edges]],{first,last}];vInLoops=Association@Map[#->EdgeCount[chain,# \[UndirectedEdge] _]&,vInLoops];vInLoops =Keys[Select[vInLoops,#==2 &]];reducedchain = VertexDelete[chain,vInLoops];If[AcyclicGraphQ[reducedchain],
Return[reducedchain]];

edges= Flatten[FindCycle[reducedchain,Infinity,All],2];res=Association@Map[#->EdgeDelete[reducedchain,#]&,edges];res=Select[res,gIsLinear];res=Select[res,gIsChain];If[Length[res]==0,Return[Missing["Can't decyclic"]]];
(*If[Length[res]>1,w=res;Print["Multi decyclics",chain,res]];*)
res=First[res];
If[AcyclicGraphQ[res],Return[res]];

Print["Can't decycle",chain];
Return[Missing[]];

];

edgePairInLROrder[chain_,a_ \[UndirectedEdge] b_] := Module[{xz},
xz = AnnotationValue[{chain,{a,b}},VertexCoordinates];
If[First[xz[[1]]]<= First[xz[[2]]],{a ,b}, {b,a}]
];

partnerNodesSortedByHeight[chain_,node_,partnerList_] := Module[{res},
res= SortBy[partnerList,Last@AnnotationValue[{chain,#},VertexCoordinates]&];
res = Map[# \[UndirectedEdge] node &,res];
res
];

deletableLeftEdges[chain_,node_] := Module[{edgeList},
edgeList= EdgeList[chain,node \[UndirectedEdge] _];
edgeList = Map[edgePairInLROrder[chain,#]&,edgeList];
edgeList = Select[edgeList,Last[#]==node&];
edgeList = Complement[Flatten[edgeList],{node}];
edgeList = partnerNodesSortedByHeight[chain,node,edgeList];
edgeList= If[Length[edgeList]>0,Drop[edgeList,-1],{}];
edgeList
];
deletableLeftEdges[chain_] := Flatten@Map[deletableLeftEdges[chain,#]&,VertexList[chain]]

deCycleChain[chain_] := Module[{reducedchain,edges},
edges= deletableLeftEdges[chain];
reducedchain= EdgeDelete[chain,edges];

If[AcyclicGraphQ[reducedchain] && chainEncircles[chain],
Return[reducedchain]];

Print["Calling old method for chain ", chain];
oldDeCycleChain[chain]

];

(*
cyclicChainParastichy[chain_] := Module[{reducedchain,v,first,last,edges,res},v=VertexList[chain];v=SortBy[v,AnnotationValue[{chain,#},VertexCoordinates]&];first=First[v];last=Last[v];
edges= Flatten[FindCycle[chain,Infinity,All],2];vInLoops=Complement[DeleteDuplicates[Flatten[List@@@edges]],{first,last}];vInLoops=Association@Map[#->EdgeCount[chain,# \[UndirectedEdge] _]&,vInLoops];vInLoops =Keys[Select[vInLoops,#==2 &]];reducedchain = VertexDelete[chain,vInLoops];If[AcyclicGraphQ[reducedchain],
Return[acyclicChainParastichy[reducedchain]]];

edges= Flatten[FindCycle[reducedchain,Infinity,All],2];res=Association@Map[#->EdgeDelete[reducedchain,#]&,edges];res=Select[res,gIsLinear];res=Select[res,gIsChain];If[Length[res]==0,Return[Missing["Can't decyclic"]]];
(*If[Length[res]>1,w=res;Print["Multi decyclics",chain,res]];
*)
res=First[res];
If[!AcyclicGraphQ[res],
Print["Can't declyce",chain];Return[Missing[]]];
acyclicChainParastichy[res]
];*)

gIsLinear[g_] := Module[{res},
res=Map[EdgeCount[g,#\[UndirectedEdge] _]&,VertexList[g]];
And@@Map[#<=2&,res]
];
gIsChain[g_] := Module[{v,first,last,graph},
v=VertexList[g];
v=SortBy[v,AnnotationValue[{g,#},VertexCoordinates]&];
first=First[v];last=Last[v];
graph=VertexDelete[g,{first,last}];
ConnectedGraphQ[graph]
];


acyclicChainParastichy[chain_]:=Module[{res},
	If[!AcyclicGraphQ[chain],Print["aCP"];Abort[]];
	res=KeySort@Counts[Map[AnnotationValue[{chain,#},EdgeStyle]&,EdgeList[chain]]];
	(* in case left or right is missing *)
	res = Append[<|RGBColor[0, 0, 1]->0,RGBColor[1, 0, 0]->0|>,res];
	res
];



(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
(*computeNodeLowerAngles[run_] := Module[{v,edges,edgePairAngle},
v=Select[VertexList[run["ContactGraph"]],bareNumberQ];
edges=Select[Association@Map[#->lowerEdges[run["ContactGraph"], #]&,v],Length[#]==2&];
edges= Association@KeyValueMap[#1->Values@KeyTake[vectorsFromNode[run,run["ContactGraph"],#1],#2]&,edges];
edgePairAngle[{e1_,e2_}] := ArcCos[(e1 . e2)/(Norm[e1]*Norm[e2])];
edges= Map[edgePairAngle,edges];
edges
]*)


(* ::Section:: *)
(*Graphics  helpers*)


(* ::Input::Initialization:: *)
(*graphToContactLines[g_,run_,leftRightColours_] := Module[{dlines,fsort,res},
dlines = Line/@Map[diskXZ[getDiskFromRun[run,#]]&,List@@@EdgeList[g],{2}];fsort[Line[{p1_,p2_}]] := If[First[p1]<First[p2],Line[{p1,p2}],Line[{p2,p1}]];dlines=Map[fsort,dlines];res = lineCylinderIntersectionColoured[#,leftRightColours]& /@dlines;
res
];

graphToContactLines[run_,leftRightColours_] := graphToContactLines[run["ContactGraph"],run,leftRightColours];

lineCylinderIntersectionColoured[Line[{{x1_,z1_},{x2_,z2_}}],leftRightColours_] :=
 Module[{slope,col},
slope= (z2-z1)/(x2-x1);
If[slope>0,col=leftRightColours["Left"],col=leftRightColours["Right"]];If[x1<-1/2,Return[
{col,
Line[{{-1/2, z2 - slope * (x2-(-1/2))},{x2,z2}}], Line[{{x1+1,z1},{1/2,z1+slope *( 1/2-(1+x1))}}]
}
]];
If[x2>1/2,Return[
{col,
Line[{{-1/2, z2- slope * (x2-1-(-1/2))},{x2-1,z2}}], Line[{{x1,z1},{1/2,z1+slope *( 1/2-(x1))}}]
}
]];
Return[{col,Line[{{x1,z1},{x2,z2}}]}] 
];
*)
diskAndVisibleCopies[Disk[{x_,z_},r_]] := {Disk[{x,z},r],If[x+r>1/2,moveDiskLeft[Disk[{x,z},r]],Nothing[]],If[x-r<-1/2,moveDiskRight[Disk[{x,z},r]],Nothing[]]};




(* ::Section:: *)
(*Extract  non - overlapping  chains*)


(* ::Input::Initialization:: *)
postRunExtractNonOverlappingChains[run_] := 
Module[{res,ppair,chainParastichies,oChains},
res=run;
oChains= nonOverlappingChains[res["RunChains"]];
ppair[redblue_] := Sort[{redblue[Red],redblue[Blue]}];
oChains =Map[Append[#,"ParastichyPair"-> ppair[#["Parastichy"]]]&,oChains];
chainParastichies= Map[Association@countsByLength[#Chain]&,oChains];
res = KeyTake[res,"Arena"];
res= Append[res,<|"ChainParastichies"->chainParastichies,"ParastichyPairs"->Map[#ParastichyPair&,oChains]|>];
res
];
nonOverlappingChains[chains_] := Module[{nextChain,chainedSoFar,thisChain,res,index},
nextChain[chain_] := Module[{w},
w=chain;
w["ChainedSoFar"]=chainedSoFar;
w["ThisChain"]=thisChain;
w["ThisChain"]=Complement[w["ThisChain"],w["ChainedSoFar"]];
w["ThisChain"] =Union[chainDiskNumbers[chain["Chain"]],w["ThisChain"]];
w["New"] =!IntersectingQ[w["ChainedSoFar"],w["ThisChain"]];
If[w["New"],
w["ChainedSoFar"]=w["ThisChain"];
w["ThisChain"]={}
];
chainedSoFar= w["ChainedSoFar"];
thisChain=w["ThisChain"];
w
];

chainedSoFar = {};thisChain={};
 res=Map[nextChain,chains];
res=Select[res,#New&];
index=1;
res= Map[Append[#,"ChainNumber"->index++]&,res];
res=Map[KeyTake[#,{"Parastichy","Chain","Radius"}]&,res];
res
];

chainDiskNumbers[chain_] := Select[VertexList[chain],bareNumberQ];



(* ::Input::Initialization:: *)
edgeDirectedRightwards[g_,v1_ \[UndirectedEdge] v2_] := Module[{v1xy,v2xy},
v1xy=AnnotationValue[{g,v1},VertexCoordinates];
v2xy=AnnotationValue[{g,v2},VertexCoordinates];
If[First[v1xy]<First[v2xy],
v1\[DirectedEdge] v2,v2\[DirectedEdge] v1]
];

directedGraphRightwards[g_] := Module[{v,vxy,edges,d,edgestyles},
v=VertexList[g];
vxy=Map[#->AnnotationValue[{g,#},VertexCoordinates]&,v];
edges =EdgeList[g];
dedges= Map[edgeDirectedRightwards[g,#]&,edges];
edgeStyles= Table[ dedges[[i]]-> AnnotationValue[{g,edges[[i]]},EdgeStyle],{i,Length[edges]}];
d=Graph[v,dedges,VertexCoordinates->vxy,VertexSize->Tiny,VertexLabels->"Name"
,EdgeShapeFunction->{{"HalfFilledArrow","ArrowSize"->.03}}
,EdgeStyle-> edgeStyles
];
d
];
strictSubgraph[g_,edges_] := Module[{res},
res=Subgraph[g,edges];
res=EdgeDelete[res,EdgeList[res]];
res= EdgeAdd[res,edges];
edgeStyles= Table[ EdgeList[res][[i]]-> AnnotationValue[{g,EdgeList[res][[i]]},EdgeStyle],{i,Length[EdgeList[res]]}];

AnnotationValue[res,EdgeStyle]=edgeStyles;
res

];

chainsByLength[chain_] := Module[{res},
dchain = directedGraphRightwards[chain] ;
nodes=SortBy[VertexList[dchain],AnnotationValue[{dchain,#},VertexCoordinates][[1]]&];

lastPath=FindShortestPath[dchain,First[nodes],Last[nodes]];
res={lastPath};
For[len=Length[lastPath]+1,len<1000,len++,
nextPath=First@FindPath[dchain,First[nodes],Last[nodes],len];
If[Length[nextPath]==Length[lastPath],Break[]];
lastPath=nextPath;
res=Append[res,lastPath];
];
res 
];
lineParastichyCount[chain_] := Module[{res},
(* chain is now aycclic *)
res=Counts[Map[AnnotationValue[{chain,#},EdgeStyle]&,EdgeList[chain]]];
res=KeySort[res];
res
];
countsByLength[chain_] := Module[{},
chains=chainsByLength[chain];
edges= Map[Apply[UndirectedEdge,#]&,Map[Partition[#,2,1]&,chains],{2}];
lines = Map[strictSubgraph[chain,#]&,edges];
lineCounts=Map[lineParastichyCount,lines];
summarize[col_] := Module[{res=Map[#[col]&,lineCounts]},
Counts[res]/Length[res]];
Map[#->summarize[#]&,{RGBColor[0, 0, 1],RGBColor[1, 0, 0]}]

];


(* ::Section:: *)
(*End*)


End[];


EndPackage[];
