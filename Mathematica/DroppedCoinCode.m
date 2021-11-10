(* ::Package:: *)

(* ::Section:: *)
(*Dropped coin model setup*)


(* ::Section:: *)
(*newdroppedcoin*)


(* ::Input::Initialization:: *)


Clear[coinD,coinZ];
coinNumber[{n_,coin_Disk}] := n;
bareCoinNumber[{n_,coin_Disk}] := bareNumber[n];
isPrincipal[coin_] := IntegerQ[coinNumber[coin]];
bareNumber[n_] := n;
bareNumber[left[n_]] := n;
bareNumber[right[n_]] := n;

lastCoinNumber[coinChain_] := Max[Select[Map[coinNumber,coinChain],IntegerQ]];

subsetCoinCollection[coinCollection_,coinNumbers_] := Select[coinCollection,MemberQ[coinNumbers,coinNumber[#] ]&];

getCoinByNumber[n_,coinCollection_] := Module[{coinSet},
coinSet = Select[coinCollection,coinNumber[#]==n &];
If[Length[coinSet]==0,Print["Can't find coin ", n];Return[{}]];
Return[First[coinSet]];
];

getCoinAndCopiesByNumber[ n_,coinCollection_] := Module[{nbare=bareNumber[n],coinset},
Select[coinCollection,MemberQ[{nbare,left[nbare],right[nbare]},coinNumber[#] ] &]
];

coinTranslateL[Disk[xy_,r_],phi_,cylinderCircumference_] := If[phi==0, Disk[xy-{cylinderCircumference,0},r],Disk[RotationTransform[2phi][xy],r]];
coinTranslateR[Disk[xy_,r_],phi_,cylinderCircumference_] :=  If[phi==0,Disk[xy+{cylinderCircumference,0},r],Disk[RotationTransform[-2phi][xy],r]];
coinTranslateL[{n_,coin_Disk},phi_,cylinderCircumference_] := {left[n],coinTranslateL[coin,phi,cylinderCircumference]};
coinTranslateR[{n_,coin_Disk},phi_,cylinderCircumference_] := {right[n],coinTranslateR[coin,phi,cylinderCircumference]};



coinH[Disk[xy_,r_]] := xy[[2]];
coinH[{n_,coin_Disk}] := coinH[coin];
coinD[Disk[xy_,r_]] := xy[[1]];
coinD[{n_,coin_Disk}] := coinD[coin];
coinR[Disk[xy_,r_]] := r;
coinR[{n_,coin_Disk}] := coinR[coin];
coinXY[{n_,Disk[xy_,r_]}] := xy;
coinZ[Disk[xy_,r_],phi_] := If[phi==0,xy[[2]],Sqrt[xy . xy]];
coinZ[{n_,coin_Disk},phi_] := coinZ[coin,phi];
coinTheta[Disk[{x_,y_},r_],phi_] := If[phi==0,2\[Pi] coinD[Disk[{x,y},r]],ArcTan[x,y]];
coinTheta[{n_,coin_Disk},phi_] := coinTheta[coin,phi];
coinDisk[{n_,coin_Disk}] := coin;
extendRadius[Disk[xy_,r_],inhibition_] := Disk[xy,r+inhibition];
extendRadius[{n_,coin_Disk},inhibition_] := extendRadius[coin,inhibition];


leftOf[n_,m_, coinCollection_,phi_] := If[phi>0,coinTheta[getCoinByNumber[n,coinCollection],phi] >  coinTheta[getCoinByNumber[m,coinCollection],phi],
 coinXY[getCoinByNumber[n,coinCollection]][[1]] < coinXY[getCoinByNumber[m,coinCollection]][[1]] ];

coinsOverZ[coinChain_,z_,phi_] :=  coinChain[[ Flatten@Position[Map[coinZ[#,phi]>z &,coinChain],True] ]];

coinDistance[Disk[xy1_,r1_],Disk[xy2_,r2_]] := Module[{v},v=xy1-xy2; Sqrt[v . v]];
coinDistance[{n1_,coin1_Disk},{n2_,coin2_Disk}] :=coinDistance[coin1,coin2];




iBoundary[coinChain_,distance_] := Map[extendRadius[#,distance]&,coinChain];


whichMinConeRegionPoint[region_,phi_] := Module[{pts,hpts,x,y,xy},
pts = MeshPrimitives[region,0];
If[phi==0,
hpts = pts /. Point[{x_,y_}] -> y,
hpts= pts /. Point[xy_] -> xy . xy
];
pts[[Ordering[hpts][[1]]]]
];


nearestCoinPair[coin_,coinChain_] := Module[{coinDistances,nearestCoins},
coinDistances = Map[coinDistance[coin,#]&,coinChain];
nearestCoins = Take[coinChain[[ Ordering[coinDistances] ]],2];
nearestCoins
];

coinLowerNeighbours[coin_,coinChain_] := Module[{nearestCoins ,res},
nearestCoins = nearestCoinPair[coin,coinChain];
If[!isPrincipal[coin],
	nearestCoins= Select[nearestCoins, coinDistance[#,coin]< 2.1 coinR[#] &]];
nearestNumbers = Map[coinNumber,nearestCoins];
res = {coinNumber[coin],nearestNumbers};


angles = Map[coinPairAngle[coin,#]&,nearestCoins];
res ={coinNumber[coin],nearestNumbers,angles};
res
];

coinPairAngle[coinm_,coinn_] := Module[{xym,xyn,v,theta},

xym= coinXY[coinm];
xyn= coinXY[coinn];
v = xyn - xym;
theta = Apply[ArcTan,v];
theta = theta/(2\[Pi]);
If[theta<0,theta = theta+1];
theta =  1/4- theta ;
If[theta<0,theta=theta+1];
theta
];

thetac[angl_] := Module[{res},res=angl- 1/2; If[res<0,res=res+1];res];



(* ::Input::Initialization:: *)
(* not used by algo, but useful for ics and explains the cone geometry *)
toCone[{d_,h_},phi_] := Module[{},
If[phi==0,
{d,h},
(h ) * { Sin[ 2 d phi] ,Cos[ 2  d phi]}]
];



coinWithConeTranslations[nextCoin_,r1_,phi_,thisCylinderCircumference_] := Module[{res},
res = {nextCoin};
If[coinD[nextCoin]+2 r1 > thisCylinderCircumference/2,
 res = Append[res,coinTranslateL[nextCoin,phi,thisCylinderCircumference]]
];
If[coinD[nextCoin] -2 r1  < -thisCylinderCircumference/2,
 res = Append[res,coinTranslateR[nextCoin,phi,thisCylinderCircumference]]
];
res
];


placeNextCoinCone[coinChain_,coinRadius_,lastCoinZ_,phi_,cylinderLU_] := Module[{ib,iboundaryMesh,cyl,availableRegion,nextdh,d,h,nextCoin,x,y},
ib = iBoundary[coinChain,coinRadius];
iboundaryMesh = DiscretizeRegion@RegionUnion[ib];
cyl = coneRegion[phi,{lastCoinZ,cylinderLU[[2]]}];
availableRegion = RegionDifference[cyl,iboundaryMesh];
nextdh  = whichMinConeRegionPoint[availableRegion,phi]; (* returns Point *) 
nextCoin =  nextdh/. (Point[{x_,y_}] ->  {lastCoinNumber[coinChain]+1,Disk[{x,y},coinRadius]});
If[Length[nextCoin!=2],Print["Can't place next coin with chain\n",coinChain];Abort[]];
nextCoin 
]


vector[m_,n_,coinCollection_] := coinXY[getCoinByNumber[n,coinCollection]] - coinXY[getCoinByNumber[m,coinCollection]]  ;

listNumbersBySteepness[n_,ncoins_,coinCollection_] := Module[{coinVectors},
coinVectors = Map[vector[n,#,coinCollection]&,ncoins];
coinVectors = Map[First,coinVectors]; (* all the same len and upper right so steepest has smallest x *)
ncoins[[ Ordering[coinVectors] ]]
];


(* ::Input:: *)
(**)
(**)


(* ::Input::Initialization:: *)
moveNumberLeft[n_] := Which[Head[n]===left,Missing[],Head[n]=== right,n[[1]],True,left[n]];
moveNumberRight[n_] := Which[Head[n]===left,n[[1]],Head[n]=== right,Missing[],True,right[n]];


addInto[res_,n_,mangle_] := Module[{r},
r = res;
If[MissingQ[r[n]],r[n]={}];
AppendTo[r[n],mangle];
r[n]=SortBy[r[n],#[[2]]&];
r
];

symmetricAdd[angleResult_,n_,m_,angle_] := Module[{res},
(* n is the newly added coin so not left or right but m might be *)
res = angleResult;
res = addInto[res,n,{m,angle}];
res = addInto[res,m,{n,thetac[angle]}];

If[Head[m]=== right,
res = addInto[res,m[[1]],{ left[n],thetac[angle]}]
];
If[Head[m]=== left,
(*Print["Adding ", n , " connected to ",m];
*)res = addInto[res,m[[1]],{ right[n],thetac[angle]}]
];
res
];



(* ::Input::Initialization:: *)

getPhi[cylinderCircumferenceFunction_,cylinderLU_] := ArcTan[
cylinderLU[[2]] - cylinderLU[[1]],cylinderCircumferenceFunction[cylinderLU[[2]] ]- cylinderCircumferenceFunction[cylinderLU[[1]]]
];

coneRegion[phi_,{h0_,h1_}] := Module[{d1,d0},
If[phi==0,
Return[DiscretizeRegion[Rectangle[{-1/2,h0},{1/2,h1}]]]
];
{d1,d0}  = Map[DiscretizeRegion,{Disk[{0,0},h1, {\[Pi]/2 - phi, \[Pi]/2 + phi}],Disk[{0,0},h0]}];
Return[RegionDifference[d1,d0]]
];



(* ::Input::Initialization:: *)
updateNeighbourFunction[nextCoinSet_,arenaAssociation_,nodeAssociation_] := Module[
{lowerN,neighbourAssociationsres,n,i,coinCollection,angles,angleResult,nextCoin},

coinCollection= nodeAssociation[Coins];

nextCoin=First@nextCoinSet;
n = coinNumber[nextCoin];

lowerN = coinLowerNeighbours[nextCoin,coinCollection][[2]];
angles= coinLowerNeighbours[nextCoin,coinCollection][[3]];

angleResult = nodeAssociation[ContactAngle];
For[i=1,i<= Length[lowerN],i++,
angleResult = symmetricAdd[angleResult,n,lowerN[[i]],angles[[i]]]
];

neighbourAssociationsres= nodeAssociation;
neighbourAssociationsres[ContactAngle]=angleResult;
neighbourAssociationsres[Coins]=Join[nodeAssociation[Coins],nextCoinSet];
neighbourAssociationsres[LastCoinZ]=coinZ[nextCoin,arenaAssociation["phi"]];


Return[neighbourAssociationsres];
];

findChainFromNodes[nextCoinSet_,arenaAssociation_,nAres_]  := Module[{coinChain2,lastChainNumbers},

lastChainNumbers =findChain[coinNumber[First[nextCoinSet]],nAres];
If[!MemberQ[lastChainNumbers,"KeyAbsent"],
(* main way *)
coinChain2 = Map[getCoinAndCopiesByNumber[#,nAres[Coins]]&,lastChainNumbers];
coinChain2 = DeleteDuplicates@Flatten[coinChain2,1]
,
(* initially or after a restart *)
coinChain2 = Join[ nAres[Coins],nextCoinSet];
];
SortBy[coinChain2,#[[2,1,1]]&]
];
findChainNumbersFromNodes[nextCoinSet_,arenaAssociation_,nAres_]  := Module[{chain},
chain = findChainFromNodes[nextCoinSet,arenaAssociation,nAres] ;
coinNumber /@ chain
];

chainFromChainNumbers[nodeAssociation_] := Module[{res},
res = Map[getCoinByNumber[#,nodeAssociation[Coins]]&,nodeAssociation[ChainNumbers]];
If[MissingQ[res],nodeAssociation[Coins],res]
];


(* ::Input::Initialization:: *)
addNextCoinCone[arenaAssociation_,nodeAssociation_] := Module[{r,nextCoin,nextCoinSet,nAres,para2,phi},

r = arenaAssociation["rFunction"][nodeAssociation[LastCoinZ]];
phi = arenaAssociation["phi"];

nextCoin =  placeNextCoinCone[  
chainFromChainNumbers[nodeAssociation]
,r
,nodeAssociation[LastCoinZ]
,phi
, arenaAssociation["cylinderLU"]
];

nextCoinSet = coinWithConeTranslations[nextCoin,r,phi,
arenaAssociation["cylinderCircumferenceFunction"][coinH[nextCoin]]];

nAres = updateNeighbourFunction[nextCoinSet,arenaAssociation,nodeAssociation];

nAres[ChainNumbers] =findChainNumbersFromNodes[nextCoinSet,arenaAssociation,nAres];

para2 = nAres[Parastichy];
para2[coinNumber[nextCoin]] =chainParastichyCount[nAres[ChainNumbers],nAres];
nAres[Parastichy] = para2;

nAres

];




(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
stackChainFromIC[icCoins_,arenaAssociation_] := Module[
{lastCoinZ,k,phi,parastichyTop,nodeAssociation},

lastCoinZ = Max[Map[coinZ[#,arenaAssociation["phi"]]&,icCoins]];


nodeAssociation =  Association[
ContactAngle->Association[]
,Parastichy->Association[]
,ChainNumbers-> coinNumber /@ icCoins
,Coins->icCoins
,LastCoinZ-> lastCoinZ];




Monitor[
For[k=0,k< arenaAssociation["coinMax"],k++,
{tim,nodeAssociation } = Timing[ 
addNextCoinCone[arenaAssociation,nodeAssociation] 
];

parastichyTop = KeyTake[nodeAssociation[Parastichy],Last[Keys@nodeAssociation[Parastichy]]];
lastCoinZ =  nodeAssociation[LastCoinZ];
If[lastCoinZ+ 2 arenaAssociation["rFunction"][lastCoinZ]> (arenaAssociation["cylinderLU"])[[2]],Break[]];
],
 {ProgressIndicator[lastCoinZ,{0,(arenaAssociation["cylinderLU"])[[2]]}],parastichyTop,tim}
];
<|"Run"-> nodeAssociation,"Arena"->arenaAssociation|>
];




(* ::Input::Initialization:: *)

chainTransitions[chain_,neighbourAssociations_] := Module[{transitions},
pairs = Partition[chain,2,1];
transitionType[{m_,n_}] := Module[{mbare=bareNumber[m],nbare=bareNumber[n]},
If[n=="KeyAbsent",Return[ Nothing]];

angles = neighbourAssociations[ContactAngle][m];

If[MissingQ[angles],Return[Nothing]];

angles = Select[angles,First[#]==n &];
If[Length[angles]==0,Return[Nothing]];
{m,n,angles[[1,2]]}
];
Map[transitionType,pairs]
];


chainParastichyCount[chain_,neighbourAssociations_] := Module[{transitions},
transitions = chainTransitions[chain,neighbourAssociations];
lr[{m_,n_,angle_}] := If[ 0 <= angle< 0.25 || 0.75 <= angle <= 1 ,Up,Down];
ud = Map[lr,transitions];
Apply[List,Counts[ud]]
];





(* ::Input::Initialization:: *)
rotateList[list_,element_] := Module[{},
If[!MemberQ[list,element],Return[list]];
pos = First@First@Position[list,element];
Join[Drop[list,pos-1],Take[list,pos-1]]
];



(* ::Input::Initialization:: *)

findChain[firstNumber_,neighbourAssociations_] := Module[{lastNumberCoinAdded,lastNumberBare,nextInChainList,nextCoinNumber,nextCoinBare,thisChain,oldchain},

thisChain = {firstNumber};

dbg = {};printDebug = MemberQ[dbg,firstNumber];


For[k=0,k<Infinity,k++,
lastNumberCoinAdded = thisChain[[ -1 ]] ;
lastNumberBare =  bareNumber[lastNumberCoinAdded ];

(* we have come in on a neighbour connection (usually but not always from the left). We want the next outward collection clockwise from that *)

allOutward = Map[First,neighbourAssociations[ContactAngle][lastNumberBare]]; 
If[Length[thisChain]==1,
If[Length[allOutward]==0,Print["No outwarda at ", firstNumber]];
nextCoinNumber = First[allOutward],
(* Length>1:  we came from somewhere *)
previousNumberCoinAdded = thisChain[[-2]];
(* reorder clockwise from there *)
If[printDebug,Print[firstNumber,":",lastNumberCoinAdded,":",allOutward]];

If[Head[lastNumberCoinAdded]==right ,
thisChain = Append[thisChain,bareNumber[lastNumberCoinAdded]];
previousNumberCoinAdded=left[bareNumber[previousNumberCoinAdded]];
];
If[Head[lastNumberCoinAdded]==left ,
If[MemberQ[Keys[neighbourAssociations[ContactAngle]],lastNumberCoinAdded],allOutward = Map[First,neighbourAssociations[ContactAngle][lastNumberCoinAdded]]]; 
];
allOutward = rotateList[allOutward,previousNumberCoinAdded];
allOutward = Cases[allOutward,Except[previousNumberCoinAdded]];
If[Length[allOutward]== 0,
(* can't finish thisChain *)
thisChain = Append[thisChain,"KeyAbsent"];
Break[]];


nextCoinNumber = First[allOutward];
];

nextCoinBare = bareNumber[nextCoinNumber ];
If[nextCoinBare == firstNumber,
thisChain = Append[thisChain,nextCoinNumber];
If[Head[nextCoinNumber]=== right,
thisChain = Append[thisChain,nextCoinBare]];
Break[]];
If[MemberQ[thisChain,nextCoinBare],Print["Chain loop from ",firstNumber, " at ",nextCoinNumber](*;Abort[]*)];
thisChain = Append[thisChain,nextCoinNumber];
];
If[printDebug,Print[thisChain]];

thisChain
];



(* ::Input::Initialization:: *)


boundingRectangle[region_] := Apply[Rectangle,Transpose[RegionBounds[region]]];
xDilate[Rectangle[{lx_,ly_},{ux_,uy_}],s_] := Rectangle[{ s lx,ly},{ s ux,uy }]


regionToPolygon[region_] :=  MeshPrimitives[BoundaryDiscretizeRegion[region],2] (* use Polygon not Region for display options *)


(* ::Section:: *)
(*Fixed cylinder*)


(* ::Input::Initialization:: *)


redConnections[n_,nANew_] := leftrightConnections[n,nANew,right];
blueConnections[n_,nANew_] := leftrightConnections[n,nANew,left];
leftrightConnections[n_,nANew_,leftright_] := Module[{connectionAngles,connections,aboveAngles,belowAngles},
connectionAngles=nANew[ContactAngle][n];
aboveAngles = Select[connectionAngles,bareNumber[n] >  bareNumber[#[[1]]]&];
If[leftright=== right ,
aboveAngles = Select[aboveAngles, #[[2]] <= 0.5  &],
aboveAngles = Select[aboveAngles, #[[2]]  > 0.5 &]
];
belowAngles = Select[connectionAngles,bareNumber[n] <  bareNumber[#[[1]]]&];
If[leftright=== right ,
belowAngles = Select[belowAngles, #[[2]]  > 0.5  &],
belowAngles = Select[belowAngles, #[[2]]  <= 0.5 &]
];

connections = Map[{n,First[#]}&,Join[aboveAngles,belowAngles]];

connections
];


(* ::Input:: *)
(**)
