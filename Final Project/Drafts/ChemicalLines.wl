(* ::Package:: *)

(* ::Section:: *)
(*Header*)


(* ::Section:: *)
(*Set up the package context, including any imports*)


BeginPackage["ChemicalRead`"]


(* ::Subsection:: *)
(*Usage Messages*)


(* ::Text:: *)
(*Find a prefix for my package*)


listPointsSlope::usage="use me"


CLListPointsSlope::usage="blah"


(* ::Subsection:: *)
(*Begin the private context*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Unprotect any system functions and local (static) variables*)


(* ::Subsection:: *)
(*Definition of auxiliary functions and local (static) variables*)


(* ::Subsection:: *)
(*Error messages for the exported objects*)


(* ::Subsection:: *)
(*Definition of the exported functions*)


(* ::Subsubsection::Closed:: *)
(*Line Utilities*)


(*Unprotect[listPointsSlope]*)
Clear[listPointsSlope]

listPointsSlope[p_?MatrixQ]:=
Quiet[(#[[2,2]]-#[[1,2]])/(#[[2,1]]-#[[1,1]])&/@Partition[N[p],2,1],{Power::infy,Infinity::indet}]


Clear[listPointsSlope]

listPointsSlope[p_?MatrixQ]:=
Quiet[(#[[2,2]]-#[[1,2]])/(#[[2,1]]-#[[1,1]])&/@Partition[N[p],2,1],{Power::infy,Infinity::indet}]


Clear[lineSlope]

lineSlope[Line[points_?MatrixQ]]:=listPointsSlope[points]
lineSlope[Line[list_/;VectorQ[list,MatrixQ]]]:=lineSlope@*Line/@list
lineSlope[lines__:{__Line}]:=lineSlope/@lines


Clear[lineDistance]

lineDistance[Line[points_?MatrixQ]]:=Total[EuclideanDistance@@@Partition[N[points],2,1]]
lineDistance[Line[list_/;VectorQ[list,MatrixQ]]]:=lineDistance@*Line/@list
lineDistance[lines:{__Line}]:=lineDistance/@lines


(* ::Subsubsection::Closed:: *)
(*Groups of Lines*)


Clear[flattenLines]

flattenLines[Line[list_/;VectorQ[list,MatrixQ]]]:=Line/@list
flattenLines[x:{__Line}]:=flattenLines/@x
(*For already list of flat lines*)
flattenLines[x:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=x
(*flattenLines[List[x___Line]]:=flattenLines*)


(* ::Subsubsection::Closed:: *)
(*Line Detection*)


Clear[doLineDetect]
doLineDetect[image_Image,tVal_,dVal_]:=
Block[{lines},
lines=ImageLines[ColorNegate[image],tVal,dVal,Method->{"Segmented"->True}];
Return[lines]
]


Clear[pickLongLines]
pickLongLines[lines_,lCutoff_]:=Select[Flatten@flattenLines@lines,lineDistance[#]>=lCutoff&]


Clear[doLongLineDetect]

doLongLineDetect[image_Image,tVal_,dVal_,lCutoff_]:=
Block[{lines,longLines},
lines=doLineDetect[image,tVal,dVal];
longLines=pickLongLines[lines,lCutoff];

Return[longLines]
]


(* ::Subsubsection::Closed:: *)
(*Text Recognize*)


Clear[testTextRecognize]
testTextRecognize[image_,level_,options___]:=
Module[{res},
res=TextRecognize[image,level,options,"BoundingBox"];
HighlightImage[image,{"Boundary",res}]
]


(* ::Subsubsection::Closed:: *)
(*Visualization*)


Clear[visualizeLineDetect]
visualizeLineDetect[image_,tVals_,dVals_]:=With[{i=ColorNegate[image]},TableForm[Table[HighlightImage[image,{Red,ImageLines[i,t,d,Method->{"Segmented"->True}]}],{t,#1},{d,#2}],TableDirections->Row,TableHeadings->{"t = "<>ToString[#]&/@#1,"d="<>ToString[#]&/@#2}]]&[tVals,dVals]
visualizeLineDetect[image_,tVals_,dVals_]:=visualizeLineDetect[image,{tVals},{dVals}]


Clear[showLines]
showLines[image_Image,lines_:{__Line}]:=HighlightImage[image,{RandomColor[],#}&/@lines]


Clear[visualizeLongLineDetect]
visualizeLongLineDetect[image_Image,tVals_List,dVals_List,lCutoff_]:=
Block[{longLines},
TableForm[
Table[
longLines=doLongLineDetect[image,t,d,lCutoff];
HighlightImage[image,{RandomColor[],#}&/@longLines],{t,#1},{d,#2}]
,TableDirections->Row
,TableHeadings->{"t = "<>ToString[#]&/@#1,"d="<>ToString[#]&/@#2}
]
]&[tVals,dVals];
visualizeLongLineDetect[image_Image,tVal_List,dVal_?NumericQ,lCutoff_?NumericQ]:=
visualizeLongLineDetect[image,tVal,{dVal},lCutoff]
visualizeLongLineDetect[image_Image,tVal_?NumericQ,dVal_List,lCutoff_?NumericQ]:=
visualizeLongLineDetect[image,{tVal},{dVal},lCutoff]
visualizeLongLineDetect[image_Image,tVal_?NumericQ,dVal_?NumericQ,lCutoff_?NumericQ]:=
visualizeLongLineDetect[image,{tVal},{dVal},lCutoff]


(* ::Subsubsection::Closed:: *)
(*Parallel lines*)


Clear[gatherParallelLines]
(*Only works on line segments*)
(*I won't make it automatically flatten lines because that's a different function.*)

(*default tolerance*)
gatherParallelLines[x:{Repeated[Line[_?MatrixQ],{2,Infinity}]}]:=gatherParallelLines[x,0.01];

gatherParallelLines[x:{Repeated[Line[_?MatrixQ],{2,Infinity}]},tol_?NumericQ]:=
Block[{combo,grouped},
combo=Transpose[{x,Flatten@*lineSlope@x}/.Indeterminate->99999999999999];
grouped=Gather[combo,Abs[#2[[2]]-#1[[2]]]<tol&];
Return[Transpose[#][[1]]&/@grouped]
]


(* ::Subsection:: *)
(*Rules for the system functions*)


(* ::Subsection:: *)
(*Restore protection of system functions*)


(* ::Subsection:: *)
(*End the private context*)


End[]


(* ::Subsection:: *)
(*Protect exported symbols*)


(* ::Subsection:: *)
(*End the package context*)


EndPackage[]
