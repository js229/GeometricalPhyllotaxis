(* ::Package:: *)

BeginPackage["PlotStylings`"]
(* provided as a default to be overwritten *)
jStyle::usage = "Association with styling defaults";

Begin["`Private`"]

jCylinderColour  = LightGreen;
jParastichyColour = <| 1 -> Red, 2 -> Blue, 3 -> Green,
4->Black |>;


fontFamily = "Gill Sans Nova Medium";
jFont[size_] := {FontFamily-> fontFamily,FontSize-> size};

jStyle = Association[
	"CylinderColour"-> jCylinderColour,
	"FontDirective" -> jFont,
	"ParastichyColour" -> jParastichyColour
	];
	
End[]

EndPackage[]


