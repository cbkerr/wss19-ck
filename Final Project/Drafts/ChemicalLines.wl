(* ::Package:: *)

(* ::Section:: *)
(*Header*)


(* ::Text:: *)
(*TODO: update default values for some of the grouping functions: findVertices, gatherParallelLines, *)


(* ::Section:: *)
(*Set up the package context, including any imports*)


BeginPackage["ChemicalRead`"]


(* ::Subsection:: *)
(*Usage Messages*)


(* ::Text:: *)
(*Package*)


(* ::Text:: *)
(*TODO : Find a prefix for naming functions*)


ChemicalRead::usage="ChemicalRead is a package for scanning chemical structures."


(* ::Text:: *)
(*Line Operations*)


listPointsSlope::usage="";
lineSlope::usage="";
lineDistance::usage="";
doLineDetect::usage="";
pickLongLines::usage="";
doLongLineDetect::usage="";
lineIntersectQ::usage="";
lineToPointDelta::usage="";


(* ::Text:: *)
(*Visualization*)


visualizeLineDetect::usage="";
showLines::usage="";
visualizeLongLineDetect::usage="";


(* ::Text:: *)
(*Line Grouping*)


gatherParallelLines::usage="";
flattenLines::usage="";
gatherNearbyLinesByCenter::usage="";
gatherIntersectingLines::usage="";


(* ::Text:: *)
(*Vertices*)


findVertices::usage="";


(* ::Text:: *)
(*Text Recognition*)


testTextRecognize::usage="";


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


(* ::Subsubsection:: *)
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


(* ::Text:: *)
(*To use Matthew's function*)


Clear[lineToPointDelta]
lineToPointDelta[l:Line[_?MatrixQ]]:={#,l[[1,2]]-#}&[l[[1,1]]]
lineToPointDelta[l:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=lineToPointDelta/@l


(* ::Text:: *)
(*From Matthew Szudzik*)


Clear[lineIntersectQ]
lineIntersectQ[a:{Repeated[_?NumericQ,2]},b:{Repeated[_?NumericQ,2]},u:{Repeated[_?NumericQ,2]},v:{Repeated[_?NumericQ,2]}]:=Quiet[Reduce[Exists[{s,t},a+s*b==u+t*v&&0<=s<=1&&0<=t<=1]],Reduce::ratnz]
lineIntersectQ[a:Line[_?MatrixQ],b:Line[_?MatrixQ]]:=Apply[lineIntersectQ,Catenate[Map[lineToPointDelta,{a,b}]]]


(* ::Subsubsection:: *)
(*Line Grouping*)


Clear[flattenLines]

flattenLines[Line[list_/;VectorQ[list,MatrixQ]]]:=Line/@list
flattenLines[x:{__Line}]:=flattenLines/@x
(*For already list of flat lines*)
flattenLines[x:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=x
(*flattenLines[List[x___Line]]:=flattenLines*)


Clear[gatherParallelLines]
(*Only works on line segments*)
(*I won't make it automatically flatten lines because that's a different function.*)

gatherParallelLines[x:{Repeated[Line[_?MatrixQ],{1,Infinity}]},tol_?NumericQ]:=
Block[{combo,grouped},
combo=Transpose[{x,Flatten@*lineSlope@x}/.Indeterminate->99999999999999];
grouped=Gather[combo,Abs[#2[[2]]-#1[[2]]]<tol&];
Return[Transpose[#][[1]]&/@grouped]
]

(*default tolerance*)
gatherParallelLines[x:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=gatherParallelLines[x,0.01];




(* ::Text:: *)
(*Gather lines by proximity. Takes a list of lines*)


Clear[gatherNearbyLinesByCenter]
gatherNearbyLinesByCenter[lines:{Repeated[Line[_?MatrixQ],{1,Infinity}]},tol_?NumericQ]:=
Gather[lines,EuclideanDistance[Mean[First[#1]],Mean[First[#2]]]<tol&]

gatherNearbyLinesByCenter[lines:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=gatherNearbyLinesByCenter[lines,10]


(* ::Text:: *)
(*Gather lines by intersecting. *)


Clear[gatherIntersectingLines]
gatherIntersectingLines[lines:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=
Block[{s2,s3,s5,s6,s7,s8},
s2=Subsets[lines,{2}];
s3=Apply[lineIntersectQ,s2,{1}];
s5=Transpose[{s2,s3}];
s6=Select[s5,#[[2]]&];
s7=Transpose[s6][[1]];
s8=Gather[s7,IntersectingQ[Part[#1],Part[#2]]&];
Return[s8]
]


(* ::Text:: *)
(*old:*)


(*Clear[gatherIntersectingLines]
gatherIntersectingLines[lines:{Repeated[Line[_?MatrixQ],{1,Infinity}]}]:=
Block[{s2,s3,s4,s5,s6,s7,s8},
s2=Subsets[#,{2}]&/@{lines};
s3=Map[Apply[lineIntersectQ],s2,{2}];
s4=Partition[Riffle[s2,s3],2];
s5=Transpose[#]&/@s4;
s6=Select[#,#[[2]]&]&/@s5;
s7=Transpose[#][[1]]&/@s6;
s8=Gather[#,IntersectingQ[Part[#1],Part[#2]]&]&/@s7;
Return[s8]
]*)


(* ::Subsubsection:: *)
(*Find vertices*)


Clear[findVertices]

findVertices[lines:{Repeated[Line[_?MatrixQ],{2,Infinity}]},tol_?NumericQ]:=
Mean/@Select[Gather[Flatten[First/@lines,1],EuclideanDistance[#1,#2]<tol&],Length[#]>=2&]

findVertices[lines:{Repeated[Line[_?MatrixQ],{2,Infinity}]}]:=findVertices[lines,3]


(* ::Text:: *)
(*Gather lines by combination of similar slope and nearby centers along the direction of the ?*)


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
(*Text Recognize*)


Clear[testTextRecognize]
testTextRecognize[image_,level_,options___]:=
Module[{res},
res=TextRecognize[image,level,options,"BoundingBox"];
HighlightImage[image,{"Boundary",res}]
]


(* ::Subsubsection:: *)
(*Shortcut Display Images*)


Clear[showLines]
showLines[image_Image,lines_:{__Line}]:=HighlightImage[image,{RandomColor[],#}&/@lines]


(* ::Subsubsection::Closed:: *)
(*Specific Visualization to find parameters*)


Clear[visualizeLineDetect]
visualizeLineDetect[image_,tVals_,dVals_]:=
With[{i=ColorNegate[image]}
,TableForm[Table[
HighlightImage[image,{Red,ImageLines[i,t,d,Method->{"Segmented"->True}]}],{t,#1},{d,#2}],TableDirections->Row,TableHeadings->{"t = "<>ToString[#]&/@#1,"d="<>ToString[#]&/@#2}]]&[tVals,dVals]
visualizeLineDetect[image_,tVals_,dVals_]:=visualizeLineDetect[image,{tVals},{dVals}]


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
