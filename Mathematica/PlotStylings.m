(* ::Package:: *)

BeginPackage["PlotStylings`"]
(* provided as a default to be overwritten *)
jStyle::usage = "Association with styling defaults";

Begin["`Private`"]




ParastichyColour = <| 
1 -> Red, 
2 -> Blue,
3 -> Green,
4->Black |>;

SetAttributes[jExport,HoldFirst];

jExport[fig_] := Module[{figname},
figname=SymbolName[Unevaluated[fig]];
(*Export[StringJoin[figname,".jpg"],fig,ImageResolution->600];
*)Export[StringJoin[figname,".pdf"],fig,ImageResolution->600];
fig
];

jStyle = Association[
	"CylinderColour"-> LightGreen,
	"FontFamily" -> "Courier",
	"ParastichyColour" -> ParastichyColour,
	"ArrowheadSpec" -> 0.02
	];

End[]

EndPackage[]


