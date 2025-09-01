(* ::Package:: *)

BeginPackage["CombinMMBE`"]
Unprotect@@Names["CombinMMBE`*"];
ClearAll@@Names["CombinMMBE`*"];



(*In this file stress fields that are calculated using the boundary element method 
are imported as SxxFinal, SyyFinal, and TxyFinal functions.

Then fields in z coordinate system are transformed to the local zq coordinate system 
and then it is rotated to the yq local coordinate system 
(SigmaxyYq, SigmaxxYq, SigmayyYq). 

*)


SEtaEta::usage=" SEtaEta: is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function."; 
eta1::usage=" "; 



Begin["`Private`"] ;

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"];
Needs["geometry`"]


If[BE==1,
SyyFinal=<<SyyFinal.dat;
SxxFinal=<<SxxFinal.dat;
TxyFinal=<<TxyFinal.dat;]





SigmaxyYq[i_]:=Module[{f},
f=TxyFinal* Cos[-2 Theta[i]]-(SxxFinal-SyyFinal) Cos[Theta[i]] Sin[Theta[i]]];

SigmaxxYq[i_]:=Module[{f},
f=SxxFinal* Cos[-Theta[i]]^2+2*TxyFinal* Cos[Theta[i]] Sin[Theta[i]]+SyyFinal* Sin[Theta[i]]^2];

SigmayyYq[i_]:=Module[{f},
f=SyyFinal* Cos[Theta[i]]^2+SxxFinal* Sin[Theta[i]]^2-TxyFinal*Sin[2 Theta[i]]];


(*SigmaEtaEtaBE[p_]:=Block[{f,f1,f2},
f1=Re[SigmaEtaEta[p]];
Xi=1.*Zeta0[p]+I eta;
f2=SEtaEta[p];
f=f1+f2]*)




SEtaEta[i_]:=Block[{f,eta1},
x = Re[zcenter[i]]+Dq[i] Cos[eta1]Cos[Theta[i]]Cosh[Zeta0[i]]-Dq[i] Sin[eta1] Sin[Theta[i]] Sinh[Zeta0[i]];
y = Im[zcenter[i]]+(Cos[eta1] Cosh[Zeta0[i]] Dq[i] Sin[Theta[i]]+Cos[Theta[i]] Dq[i] Sin[eta1] Sinh[Zeta0[i]]);

f=(2.0/(Cosh[2*Zeta0[i]]-Cos[2*eta1]))(((Sin[eta1]*Cosh[Zeta0[i]])^2)*SigmaxxYq[i]+((Cos[eta1]*Sinh[Zeta0[i]])^2)*SigmayyYq[i]-0.5*(Sin[2*eta1]*Sinh[2*Zeta0[i]])*SigmaxyYq[i])]


End[] 
Protect@@Names["CombinMMBE`*"]
EndPackage[] 

(*BeginPackage["CombinMMBE`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","PackageManipulations`","run1`"}] 
EndPackage[] *)

(*BeginPackage["CombinMMBE`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","PackageManipulations`","run1`"}] 
EndPackage[]*)

BeginPackage["CombinMMBE`",{"geometry`"}] 
EndPackage[]











































































