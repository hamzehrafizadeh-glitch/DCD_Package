(* ::Package:: *)

BeginPackage["CombinDCDTEMP0`"]
Unprotect@@Names["CombinDCDTEMP0`*"];
ClearAll@@Names["CombinDCDTEMP0`*"];



(*In this file thermal stress fields are calculated for circular shape inclusions. 

Then fields in z coordinate system are transformed to the local zq coordinate system 
and then it is rotated to the yq local coordinate system 
(SigmaxyYq, SigmaxxYq, SigmayyYq). 

SEtaEta: is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function.*)






Begin["`Private`"] ;




(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"];
Needs["geometry`"]
Print["CombinDCDTEMP0 = ",TS];


(*The Piecewise function is more correct but it should be computationaly more expensive! That is why I am not using it.*)
SigmarrQ[i_,TStart_,TEnd_,ThermalExpansion0t_,ThermalExpansionQt_,xTemp_,yTemp_]:=Block[{radious,r,dalpha,dtem,p,f},
radious=l2[i]/(1+2*delta/l1[i]);Clear[xTemp,yTemp];

r=Sqrt[xTemp^2 + yTemp^2];
dalpha=ThermalExpansion0t-ThermalExpansionQt;
dtem=TEnd-TStart;
p=If[ThermalExpansionQt==0.0,0.0,1.*dalpha*dtem/((1+Nu0)/E0+(1-Nu[i])/EQ[i])];
Print["p =  ",p];
f=-p*(radious/r)^2.
(*Piecewise[{{p,r<radious},{f,r>=radious}}]*)
]

SigmattQ[i_,TStart_,TEnd_,ThermalExpansion0t_,ThermalExpansionQt_,xTemp_,yTemp_]:=-SigmarrQ[i,TStart,TEnd,ThermalExpansion0t,ThermalExpansionQt,xTemp,yTemp];




(*-(25/(xTemp^2+yTemp^2)^1.`)*)
If[TS==1,
SrrFinal[xTemp_,yTemp_,i_]:=-50/(xTemp^2+yTemp^2);
SrrFinalt[xTemp_,yTemp_,i_]:=SigmarrQ[i,TStart,TEnd,ThermalExpansion0,ThermalExpansionQ[i],xTemp,yTemp];

SttFinal[xTemp_,yTemp_,i_]:=-SrrFinal[xTemp,yTemp,i];
(*These are stress fields in the global coordinate system.*)



SxxFinalz0[xTemp_,yTemp_,i_]:=SrrFinal[xTemp,yTemp,i](-yTemp^2+ xTemp^2)/(xTemp^2+yTemp^2);

TxyFinalz0[xTemp_,yTemp_,i_]:=(2*SrrFinal[xTemp,yTemp,i] xTemp yTemp)/(xTemp^2+yTemp^2);

SyyFinalz0[xTemp_,yTemp_,i_]:=SrrFinal[xTemp,yTemp,i](yTemp^2- xTemp^2)/(xTemp^2+yTemp^2) ;]



SxxFinal0[xx_,yy_,i_]:=Block[{f,r,xxt,yyt},r=SxxFinalz0[xxt,yyt,i];
xxt=xx-0.Re[zcenter[i]];yyt=yy-0.Im[zcenter[i]];
r]

SyyFinal0[xx_,yy_,i_]:=Block[{f,r,xxt,yyt},r=SyyFinalz0[xxt,yyt,i];
xxt=xx-0.Re[zcenter[i]];yyt=yy-0.Im[zcenter[i]];
r]

TxyFinal0[xx_,yy_,i_]:=Block[{f,r,xxt,yyt},r=TxyFinalz0[xxt,yyt,i];
xxt=xx-0.Re[zcenter[i]];yyt=yy-0.Im[zcenter[i]];
r]

SxxFinal=Simplify[Sum[SxxFinal0[xx,yy,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]
SyyFinal=Simplify[Sum[SyyFinal0[xx,yy,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]
TxyFinal=Simplify[Sum[TxyFinal0[xx,yy,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]


(*USE THIS ONE FOR TEST, AND FOR 1 RS FIELD*)
(*rtx=0;rty=-6.;
SxxFinal=Simplify[Sum[SxxFinalz0[xx-rtx,yy-rty,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]
SyyFinal=Simplify[Sum[SyyFinalz0[xx-rtx,yy-rty,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]
TxyFinal=Simplify[Sum[TxyFinalz0[xx-rtx,yy-rty,i],{i, Times@@ Dimensions[ThermalExpansiontemp]}]]*)


SyyFinal >> SyyFinalTS.dat;
SxxFinal >> SxxFinalTS.dat;
TxyFinal >> TxyFinalTS.dat;







End[] 
Protect@@Names["CombinDCDTEMP0`*"]
EndPackage[] 

BeginPackage["CombinDCDTEMP0`",{"geometry`"}] 
EndPackage[] 


































































































