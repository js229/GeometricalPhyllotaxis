(* ::Package:: *)

BeginPackage["TextbookStylings`"]

jExpo::usage = "Exports a figure to pdf and jpg";
jFont::usage = "Typeface function";
jCylinderColour::usage  = "";
jParastichyColour::usage = "Parastichy colours";
jStyle::usage = "Association with styling defaults";

Begin["`Private`"]
jCylinderColour  = LightGreen;
jParastichyColour = <| 1 -> Red, 2 -> Blue, 3 -> Green,
4->Black |>;

SetAttributes[jExpo,HoldFirst];
jExpo[fig_] := Module[{figname},
figname=SymbolName[Unevaluated[fig]];
Export[StringJoin[(*"mma-",*)figname,".jpg"],fig,ImageResolution->600];
Export[StringJoin[(*"mma-",*)figname,".pdf"],fig,ImageResolution->600];
fig
];

fontFamily = "Gill Sans Nova Medium";
fontStyle14 = {FontFamily-> fontFamily,FontSize-> 14};
fontStyle12 = {FontFamily-> fontFamily,FontSize-> 12};
jFont[size_] := {FontFamily-> fontFamily,FontSize-> size};

jStyle = Association[
	"CylinderColour"-> jCylinderColour,
	"FontDirective" -> jFont,
	"ParastichyColour" -> jParastichyColour
	];
	
End[]

EndPackage[]

(*SetOptions[ParametricPlot,AxesStyle->Directive[FontFamily->fontFamily],
BaseStyle->Directive[fontFamily],
LabelStyle->Directive[fontFamily]];
SetOptions[Plot,AxesStyle->Directive[Red,FontFamily->fontFamily],
BaseStyle->Directive[Blue,fontFamily],
LabelStyle->Directive[Blue,fontFamily]];
*)



