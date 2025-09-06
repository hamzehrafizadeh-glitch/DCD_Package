(* ::Package:: *)

BeginPackage["CombinDCDTEMP`"]
Unprotect@@Names["CombinDCDTEMP`*"];
ClearAll@@Names["CombinDCDTEMP`*"];



(*In this file thermal stress fields are calculated for circular shape inclusions. 

Then fields in z coordinate system are transformed to the local zq coordinate system 
and then it is rotated to the yq local coordinate system 
(SigmaxyYq, SigmaxxYq, SigmayyYq). 

SEtaEta: is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function.*)



SEtaEtaTS::usage="is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function. "; 
eta2::usage=" "; 
SxxFinalz::usage=" These are stress fields in the global coordinate system."; 
TxyFinalz::usage="These are stress fields in the global coordinate system. "; 
SyyFinalz::usage="These are stress fields in the global coordinate system. "; 
SxxFinalTemp::usage=" These are stress fields in the zq coordinate system."; 
TxyFinalTemp::usage="These are stress fields in the zq coordinate system. "; 
SyyFinalTemp::usage="These are stress fields in the zq coordinate system. "; 
SigmaxyYq::usage=" These are stress fields that are rotated acording to the local ith(yi) coordinate system."; 
SigmaxxYq::usage=" These are stress fields that are rotated acording to the local ith(yi) coordinate system."; 
SigmayyYq::usage=" These are stress fields that are rotated acording to the local ith(yi) coordinate system."; 
SIFTempTemp::usage=".";





Begin["`Private`"] ;




(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"];
Needs["geometry`"]
Print["TS = ",TS];


(*(*The Piecewise function is more correct but it should be computationaly more expensive! That is why I am not using it.*)
SigmarrQ[i_,TStart_,TEnd_,ThermalExpansion0t_,ThermalExpansionQt_,xTemp_,yTemp_]:=Block[{radious,r,dalpha,dtem,p,f},
radious=l2[i]/(1+2*delta/l1[i]);Clear[xTemp,yTemp];

r=Sqrt[(xTemp-Re[zcenter[i]])^2 + (yTemp-Im[zcenter[i]])^2];
dalpha=ThermalExpansion0t-ThermalExpansionQt;
dtem=TEnd-TStart;
p=-2.*dalpha*dtem/((1+Nu0)/E0+(1-Nu[i])/EQ[i]);
(*Print["p =  ",p];*)
f=-p*(radious/r)^2.
(*Piecewise[{{p,r<radious},{f,r>=radious}}]*)
]

SigmattQ[i_,TStart_,TEnd_,ThermalExpansion0t_,ThermalExpansionQt_,xTemp_,yTemp_]:=-SigmarrQ[i,TStart,TEnd,ThermalExpansion0t,ThermalExpansionQt,xTemp,yTemp];

*)


(*(*-(516.26668151184`/(xTemp^2+yTemp^2)^1.`)*)
If[TS==1,
SrrFinal[xTemp_,yTemp_]:=Sum[SigmarrQ[i,TStart,TEnd,ThermalExpansion0,ThermalExpansionQ[i],xTemp,yTemp],{i, Times@@ Dimensions[ThermalExpansiontemp]}];
SttFinal[xTemp_,yTemp_]:=-SrrFinal[xTemp,yTemp];
(*These are stress fields in the global coordinate system.*)

SxxFinalz[xTemp_,yTemp_]:=(SttFinal[xTemp,yTemp] xTemp^2+SrrFinal[xTemp,yTemp] yTemp^2)/(xTemp^2+yTemp^2);
TxyFinalz[xTemp_,yTemp_]:=((SrrFinal[xTemp,yTemp]-SttFinal[xTemp,yTemp]) xTemp yTemp)/(xTemp^2+yTemp^2);

SyyFinalz[xTemp_,yTemp_]:=(SrrFinal[xTemp,yTemp] xTemp^2+SttFinal[xTemp,yTemp] yTemp^2)/(xTemp^2+yTemp^2);]
*)


(*cos=y/Abs[x^2 +y^2];sin=x/Abs[x^2 +y^2];
trans=cos	-sin
sin	cos

.SrrFinal	0
0	SttFinal

.cos	sin
-sin	cos

*)
If[TS==1,
fyy = <<SyyFinalTS.dat;
fxy = <<TxyFinalTS.dat;
fxx= <<SxxFinalTS.dat,
fyy = 0.;
fxy = 0.;
fxx= 0.;]





SyyFinalz[xTemp_,yTemp_]:= Module[{f},Clear[xx,yy];f = fyy;
xx=xTemp;yy=yTemp;
f];

TxyFinalz[xTemp_,yTemp_]:= Module[{f},Clear[xx,yy];f = fxy;
xx=xTemp;yy=yTemp;
f];


SxxFinalz[xTemp_,yTemp_]:= Module[{f},Clear[xx,yy];f = fxx;
xx=xTemp;yy=yTemp;
f];



(*These are stress fields in the qth zq coordinate systems. *)

SxxFinal[q_,xTemp_,yTemp_]:=Module[{f,f1,xq,yq,xt,yt},
f=SxxFinalz[xt,yt];
xt=xq+Re[zcenter[q]];
yt=yq+Im[zcenter[q]];
f1=f;
xq=xTemp;yq=yTemp;
f1]

TxyFinal[q_,xTemp_,yTemp_]:=Module[{f,f1,xq,yq,xt,yt},
f=TxyFinalz[xt,yt];
xt=xq+Re[zcenter[q]];
yt=yq+Im[zcenter[q]];
f1=f;
xq=xTemp;yq=yTemp;
f1]

SyyFinal[q_,xTemp_,yTemp_]:=Module[{f,f1,xq,yq,xt,yt},
f=SyyFinalz[xt,yt];
xt=xq+Re[zcenter[q]];
yt=yq+Im[zcenter[q]];
f1=f;
xq=xTemp;yq=yTemp;
f1]




SxxFinalTemp[q_,xTemp_,yTemp_]:=SxxFinal[q,xTemp,yTemp];
TxyFinalTemp[q_,xTemp_,yTemp_]:=TxyFinal[q,xTemp,yTemp];

SyyFinalTemp[q_,xTemp_,yTemp_]:=SyyFinal[q,xTemp,yTemp];




(*These are stress fields that are rotated acording to the local ith(yi) coordinate system.*)
SigmaxyYq[i_,xTemp_,yTemp_]:=Module[{f},
f=TxyFinal[i,xTemp,yTemp]* Cos[-2 Theta[i]]-(SxxFinal[i,xTemp,yTemp]-SyyFinal[i,xTemp,yTemp]) Cos[Theta[i]] Sin[Theta[i]]];

SigmaxxYq[i_,xTemp_,yTemp_]:=Module[{f},
f=SxxFinal[i,xTemp,yTemp]* Cos[-Theta[i]]^2+2*TxyFinal[i,xTemp,yTemp]* Cos[Theta[i]] Sin[Theta[i]]+SyyFinal[i,xTemp,yTemp]* Sin[Theta[i]]^2];

SigmayyYq[i_,xTemp_,yTemp_]:=Module[{f},
f=SyyFinal[i,xTemp,yTemp]* Cos[Theta[i]]^2+SxxFinal[i,xTemp,yTemp]* Sin[Theta[i]]^2-TxyFinal[i,xTemp,yTemp]*Sin[2 Theta[i]]];


SEtaEtaTS[i_]:=Block[{f,eta2,xTemp1,yTemp1},
f=(2.0/(Cosh[2*Zeta0[i]]-Cos[2*eta2]))(((Sin[eta2]*Cosh[Zeta0[i]])^2)*SigmaxxYq[i,xTemp1,yTemp1]+((Cos[eta2]*Sinh[Zeta0[i]])^2)*SigmayyYq[i,xTemp1,yTemp1]-0.5*(Sin[2*eta2]*Sinh[2*Zeta0[i]])*SigmaxyYq[i,xTemp1,yTemp1]);
xTemp1=Re[dd[i]Cosh[Zeta0[i]+I*eta2]];
yTemp1=Im[dd[i]Cosh[Zeta0[i]+I*eta2]];
f]
Print["TS = ",TS];


SIFTempTemp[i_]:=Block[{f,xini,xtemp,ytemp},xini=Dq[i];f=Simplify[Sqrt[2Pi*((xtemp-Re[xini])^2 +(ytemp-Im[xini])^2)](SigmayyYq[i,xtemp,ytemp]+I*SigmaxyYq[i,xtemp,ytemp])];
xtemp=Re[(EndPoint[i]-zcenter[i])Exp[-I*Thetaprim[i]]];
ytemp=Im[(EndPoint[i]-zcenter[i])Exp[-I*Thetaprim[i]]];
f]


End[] 
Protect@@Names["CombinDCDTEMP`*"]
EndPackage[] 

BeginPackage["CombinDCDTEMP`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","PackageManipulations`","run1`"}] 
EndPackage[] 
(*BeginPackage["CombinDCDTEMP`",{"geometry`"}] 
EndPackage[] *)

(*,"SIF`"*)






















































































