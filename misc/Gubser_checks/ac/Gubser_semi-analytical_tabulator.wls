#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
(* Gubser semi-analytical solution tabulator - DOI:10.1103/PhysRevC.91.014903 *)


(* ::Input::Initialization:: *)
Print["Gubser semi-analytical solution tabulator - DOI:10.1103/PhysRevC.91.014903
University of Illinois at Urbana-Champaign, 2023
"]
Print[">> This program tabulates the solutions for temperature, velocity and shear-stress tensor on a x-y spatial grid for Gubser flow with constant shear viscosity for a given hyperbolic time \[Tau]."]


(* ::Input:: *)
(*(* Instructions on how to use the script *)*)


(* ::Input::Initialization:: *)
Clear[\[Tau],gridMax,gridStep,\[HBar]c,\[Eta]OverS,\[Alpha],c,sols,\[Rho],T]
args=Rest[$CommandLine];
If[Length[args]!=5,
Print["
	Usage: mathematica/math/MathKernel -script Gubser_semi-analytical_tabulator.wls <tau [fm/c]> <grid max [fm]> <grid step [fm]>"];
Quit[]
];
\[Tau]=ToExpression[args[[3]]];
gridMax=ToExpression[args[[4]]];
gridStep=ToExpression[args[[5]]];


(* ::Input::Initialization:: *)
(* Constants *)


(* ::Input::Initialization:: *)
\[HBar]c=(197.326980/1000); 
\[Eta]OverS=0.2;
\[Alpha]=3 \[Pi]^2/90 (2(Nc^2-1)+7/2 Nc Nf)/.{Nc->3,Nf->2.5};
c=5;


(* ::Input::Initialization:: *)
(* Differential equations *)


(* ::Input::Initialization:: *)
sols=NDSolve[{
T\[LetterSpace]hat'[\[Rho]]/T\[LetterSpace]hat[\[Rho]]+2/3 Tanh[\[Rho]]==1/3 \[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta][\[Rho]]Tanh[\[Rho]],
c/T\[LetterSpace]hat[\[Rho]] \[Eta]OverS(\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta]'[\[Rho]]+4/3 (\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta][\[Rho]])^2 Tanh[\[Rho]])+\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta][\[Rho]]==4/3 \[Eta]OverS/T\[LetterSpace]hat[\[Rho]] Tanh[\[Rho]],
T\[LetterSpace]hat[0]==1.2,
\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta][0]==0},
{T\[LetterSpace]hat,\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta]},
{\[Rho],-12,12}
];


(* ::Input::Initialization:: *)
(* Change of coordinates *)


(* ::Input::Initialization:: *)
\[Rho][\[Tau]_,r_]:=ArcSinh[-((1-\[Tau]^2+r^2)/(2\[Tau]))];
T[\[Tau]_,r_]:=\[HBar]c ((T\[LetterSpace]hat[\[Rho][\[Tau],r]]/.sols)/\[Tau])[[1]];
T[\[Tau]_,x_,y_]:= T[\[Tau],Sqrt[x^2+y^2]];
u\[LetterSpace]r[\[Tau]_,r_]:=Sinh[ArcTanh[(2\[Tau] r)/(1+\[Tau]^2+r^2)]];
u\[LetterSpace]x[\[Tau]_,x_,y_]:=Piecewise[{{x/Sqrt[x^2+y^2] Sinh[ArcTanh[(2\[Tau] Sqrt[x^2+y^2])/(1+\[Tau]^2+x^2+y^2)]],Sqrt[x^2+y^2]!=0},{0.0,Sqrt[x^2+y^2]==0}}];
u\[LetterSpace]y[\[Tau]_,x_,y_]:=Piecewise[{{y/Sqrt[x^2+y^2] Sinh[ArcTanh[(2\[Tau] Sqrt[x^2+y^2])/(1+\[Tau]^2+x^2+y^2)]],Sqrt[x^2+y^2]!=0},{0.0,Sqrt[x^2+y^2]==0}}];
\[Pi]\[LetterSpace]hat\[LetterSpace]\[Eta]\[Eta][\[Tau]_,r_]:=4/3 \[Alpha] \[HBar]c ((T\[LetterSpace]hat[\[Rho][\[Tau],r]]/.sols)[[1]])^4 (\[Pi]\[LetterSpace]bar\[LetterSpace]\[Eta]\[Eta][\[Rho][\[Tau],r]]/.sols)[[1]]; 
\[Pi]\[LetterSpace]\[Eta]\[Eta][\[Tau]_,x_,y_]:=1/\[Tau]^6 \[Pi]\[LetterSpace]hat\[LetterSpace]\[Eta]\[Eta][\[Tau],Sqrt[x^2+y^2]];
\[Pi]\[LetterSpace]xx[\[Tau]_,x_,y_]:=-(1/2) 1/\[Tau]^4 (1+Piecewise[{{(x/Sqrt[x^2+y^2])^2 Sinh[ArcTanh[(2\[Tau] Sqrt[x^2+y^2])/(1+\[Tau]^2+x^2+y^2)]]^2,x^2+y^2!=0},{0.0,x^2+y^2==0}}])\[Pi]\[LetterSpace]hat\[LetterSpace]\[Eta]\[Eta][\[Tau],Sqrt[x^2+y^2]];
\[Pi]\[LetterSpace]yy[\[Tau]_,x_,y_]:=-(1/2) 1/\[Tau]^4 (1+Piecewise[{{(y/Sqrt[x^2+y^2])^2 Sinh[ArcTanh[(2\[Tau] Sqrt[x^2+y^2])/(1+\[Tau]^2+x^2+y^2)]]^2,x^2+y^2!=0},{0.0,x^2+y^2==0}}])\[Pi]\[LetterSpace]hat\[LetterSpace]\[Eta]\[Eta][\[Tau],Sqrt[x^2+y^2]];
\[Pi]\[LetterSpace]xy[\[Tau]_,x_,y_]:=Piecewise[{{-(1/2) 1/\[Tau]^4 (x y )/(x^2+y^2) Sinh[ArcTanh[(2\[Tau] Sqrt[x^2+y^2])/(1+\[Tau]^2+x^2+y^2)]]^2 \[Pi]\[LetterSpace]hat\[LetterSpace]\[Eta]\[Eta][\[Tau],Sqrt[x^2+y^2]],x^2+y^2!=0},{0.0,x^2+y^2==0}}];



(* ::Input::Initialization:: *)
(* Create output files *)


(* ::Input::Initialization:: *)
(* Table with columns "x [fm], y [fm], T [GeV], u^x, u^y, \[Pi]^xx [GeV/fm^3], \[Pi]^yy [GeV/fm^3], \[Pi]^xy [GeV/fm^3], \[Tau]^2\[Pi]^\[Eta]\[Eta] [GeV/fm^3]" as a function of time \[Tau] *)
formatted\[LetterSpace]Profile\[LetterSpace]Table[\[Tau]_:1,gridMax_:5,gridStep_:0.05,path_:Directory[]]:=
Module[
{xMax=Abs@gridMax,
xMin=-Abs@gridMax,
yMax=Abs@gridMax,
yMin=-Abs@gridMax,
\[CapitalDelta]x=Abs@gridStep,
\[CapitalDelta]y=Abs@gridStep},
outputTable=ArrayFlatten[
Table[{
DecimalForm@N[x,7],
DecimalForm@N[y,7],
DecimalForm@N[T[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]x[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]y[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xx[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]yy[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xy[\[CapitalTau],x,y],7],
DecimalForm@N[(\[CapitalTau]^2 \[Pi]\[LetterSpace]\[Eta]\[Eta][\[CapitalTau],x,y]),7]
}/.{\[CapitalTau]->\[Tau]},
{x,xMin,xMax,\[CapitalDelta]x},{y,yMin,yMax,\[CapitalDelta]y}],1];
format[x_]:=ToString@NumberForm[x,{7,6},NumberFormat->(StringJoin[{#1,""}/. "":>""/;#2==""]&)];
formattedOutput=Map[format,outputTable,{2}];
filename="Profile_tau="<>ToString[\[Tau]]<>"fm.dat";
completePath=FileNameJoin[{path,filename}];
Export[completePath,formattedOutput,"Table","FieldSeparators"->"     "];
];
formatted\[LetterSpace]yeq0\[LetterSpace]Table[\[Tau]_:1,gridMax_:5,gridStep_:0.05,path_:Directory[]]:=
Module[
{xMax=Abs@gridMax,
xMin=-Abs@gridMax,
\[CapitalDelta]x=Abs@gridStep},
outputTable=ArrayFlatten[
Table[{
DecimalForm@N[x,7],
DecimalForm@N[y,7],
DecimalForm@N[T[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]x[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]y[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xx[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]yy[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xy[\[CapitalTau],x,y],7],
DecimalForm@N[(\[CapitalTau]^2 \[Pi]\[LetterSpace]\[Eta]\[Eta][\[CapitalTau],x,y]),7]
}/.{\[CapitalTau]->\[Tau],y->0},
{x,xMin,xMax,\[CapitalDelta]x}],0];
format[x_]:=ToString@NumberForm[x,{7,6},NumberFormat->(StringJoin[{#1,""}/. "":>""/;#2==""]&)];
formattedOutput=Map[format,outputTable,{2}];
filename="y=0_tau="<>ToString[\[Tau]]<>"_SemiAnalytic.dat";
completePath=FileNameJoin[{path,filename}];
Export[completePath,formattedOutput,"Table","FieldSeparators"->"     "];
];
formatted\[LetterSpace]yeqx\[LetterSpace]Table[\[Tau]_:1,gridMax_:5,gridStep_:0.05,path_:Directory[]]:=
Module[
{xMax=Abs@gridMax,
xMin=-Abs@gridMax,
\[CapitalDelta]x=Abs@gridStep},
outputTable=ArrayFlatten[
Table[{
DecimalForm@N[x,7],
DecimalForm@N[y,7],
DecimalForm@N[T[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]x[\[CapitalTau],x,y],7],
DecimalForm@N[u\[LetterSpace]y[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xx[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]yy[\[CapitalTau],x,y],7],
DecimalForm@N[\[Pi]\[LetterSpace]xy[\[CapitalTau],x,y],7],
DecimalForm@N[(\[CapitalTau]^2 \[Pi]\[LetterSpace]\[Eta]\[Eta][\[CapitalTau],x,y]),7]
}/.{\[CapitalTau]->\[Tau],y->x},
{x,xMin,xMax,\[CapitalDelta]x}],0];
format[x_]:=ToString@NumberForm[x,{7,6},NumberFormat->(StringJoin[{#1,""}/. "":>""/;#2==""]&)];
formattedOutput=Map[format,outputTable,{2}];
filename="y=x_tau="<>ToString[\[Tau]]<>"_SemiAnalytic.dat";
completePath=FileNameJoin[{path,filename}];
Export[completePath,formattedOutput,"Table","FieldSeparators"->"     "];
]


(* ::Input::Initialization:: *)
(* Write output files *)


(* ::Input::Initialization:: *)
Print[">> Table format is...
	x [fm], y [fm], T [GeV], u^x, u^y,
	\[Pi]^{xx} [GeV/fm^3], \[Pi]^{yy} [GeV/fm^3], \[Pi]^{xy} [GeV/fm^3], \[Tau]^2 \[Pi]^{\[Eta]\[Eta]} [GeV/fm^3]"]
Print[">> Generating profile..."]
formatted\[LetterSpace]Profile\[LetterSpace]Table[\[Tau],gridMax,gridStep]
Print[">> Complete."]


(* ::Input::Initialization:: *)
Print[">> Generating slice y=0 fm..."]
formatted\[LetterSpace]yeq0\[LetterSpace]Table[\[Tau],gridMax,gridStep]
Print[">> Complete."]


(* ::Input::Initialization:: *)
Print[">> Generating slice y=x..."]
formatted\[LetterSpace]yeqx\[LetterSpace]Table[\[Tau],gridMax,gridStep]
Print[">> Complete."]
