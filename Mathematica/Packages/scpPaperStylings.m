(* ::Package:: *)

BeginPackage["TextbookStylings`"]


jExport::usage = "Exports a named object to a file with the same name";
jStyle::usage = "Association with styling defaults";



(* ::Input::Initialization:: *)
Begin["`Private`"]


CylinderColour =RGBColor[0.8692074837412045, 0.86, 0.6466894826847229];
CylinderColour3D =Blend[{CylinderColour,Green},0.2];
ParastichyColour = <|1->RGBColor[0.7194116276462343, 0.32145384723880266`, 0.27090344303450226`],2->RGBColor[0.4687348547943994, 0.29473626986998797`, 0.5955133155492501], 3-> RGBColor[0.4518714536621338, 0.549629484806721, 0.24743503166788056`], 4->RGBColor[0.955983350983375, 0.9459627980758899, 0.942555568392974],5->RGBColor[0.2138260590152513, 0.18751247299324145`, 0.15071362998820095`]|>;


SetAttributes[jExport,HoldFirst];

pageWidthPT  = 341; (* from latex layout *) 

jExport[fig_,size_:1] := Module[{figname,imageWidth},

figname=SymbolName[Unevaluated[fig]];
imageWidth = pageWidthPT * size;
Export[StringJoin[figname,".jpg"],fig,ImageWidth->imageWidth,ImageResolution->300];
Export[StringJoin[figname,".pdf"],fig,ImageWidth->imageWidth,ImageResolution->300];
fig
];




jStyle = Association[
	"CylinderColour"-> CylinderColour,
	"CylinderColour3D"-> CylinderColour3D,
	"FontFamily" -> "Gill Sans MT",
	"FullWidthImagePT" -> pageWidthPT,
	"LaTeXBodyFamily" -> "EB Garamond 12", (* but gives old-style numbers *)
	"LaTeXDefaultFamily" -> "Latin Modern Roman 12", (* gives inline numbers *)
	"ParastichyColour" -> ParastichyColour,
	"ArrowheadSpec" -> 0.02
	];

End[]

EndPackage[]








