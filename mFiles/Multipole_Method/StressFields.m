(* ::Package:: *)

BeginPackage["StressFields`"]
Unprotect@@Names["StressFields`*"];
ClearAll@@Names["StressFields`*"];


SigmaEtaEta::usage=" \[Sigma]\[Eta]\[Eta]; stress field in yq local coordinate system"; 
Sigma11m::usage=" sigma11; stress field in yq local coordinate system"; 
Sigma12m::usage=" sigma12; stress field in yq local coordinate system"; 
Sigma22m::usage=" sigma22; stress field in yq local coordinate system"; 
Sigma11z::usage=" sigma11; stress field in zq local coordinate system"; 
Sigma12z::usage=" sigma12; stress field in zq local coordinate system"; 
Sigma22z::usage=" sigma22; stress field in zq local coordinate system"; 
Sigma1111z::usage=" sigma11; stress field in z global coordinate system"; 
Sigma1211z::usage=" sigma12; stress field in z global coordinate system"; 
Sigma2211z::usage=" sigma22; stress field in z global coordinate system"; 
PlotSigma11zQ::usage=" plot sigma11; stress field in zq local coordinate system"; 
PlotSigma12zQ::usage=" plot sigma12; stress field in zq local coordinate system"; 
PlotSigma22zQ::usage=" plot sigma22; stress field in zq local coordinate system"; 
PlotSigma11z::usage=" plot sigma11; stress field in z global coordinate system"; 
PlotSigma12z::usage=" plot sigma11; stress field in z global coordinate system"; 
PlotSigma22z::usage=" plot sigma11; stress field in z global coordinate system"; 
x::usage=" "; 
y::usage=" "; 


Sigma11mStress::usage=" sigma11; stress field in yq local coordinate system"; 
Sigma12mStress::usage=" sigma12; stress field in yq local coordinate system"; 
Sigma22mStress::usage=" sigma22; stress field in yq local coordinate system"; 
Sigma11zStress::usage=" sigma11; stress field in zq local coordinate system"; 
Sigma12zStress::usage=" sigma12; stress field in zq local coordinate system"; 
Sigma22zStress::usage=" sigma22; stress field in zq local coordinate system"; 
Sigma1111zStress::usage=" sigma11; stress field in z global coordinate system"; 
Sigma1211zStress::usage=" sigma12; stress field in z global coordinate system"; 
Sigma2211zStress::usage=" sigma22; stress field in z global coordinate system"; 


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
(* ***Get["PhiPsi.m"];*)
Get["PhiPsi.m"];
Needs["PhiPsi`"];
Needs["CombinDCDTEMP`"];

Print["I am in Reloading files"]


(*SigmaEtaEta[q_], Sigma11mStress[q_], Sigma12mStress[q_], Sigma22mStress[q_], and 

Sigma11mStresslist,
Sigma12mStresslist,
Sigma22mStresslist,

are stress terms due to interactions. *)

SigmaEtaEta[q_]:=Block[{f},f=(Beta0[q]+BetaFar[q])-Re[Conjugate[Tau0[q]+TauFar[q]]]]
Unprotect[z,Xi,x,y];
(*z = x + I y*)
Sigma11mStress[q_]:=Module[{f},f=Re[AlphaFar[q]+Alpha0[q]]]
Sigma12mStress[q_]:=Module[{f},f=-Im[AlphaFar[q]+Alpha0[q]]]

Sigma11mStresslist:=Module[{f},f=Table[Sigma11mStress[q],{q,ntot}]]
Sigma12mStresslist:=Module[{f},f=Table[Sigma12mStress[q],{q,ntot}]]

Sigma22mStress[q_]:=Module[{f},f=BetaFar[q]+Beta0[q]-Re[Alpha0[q]]-Re[AlphaFar[q]]]
Sigma22mStresslist:=Module[{f},f=Table[Sigma22mStress[q],{q,ntot}]]



Sigma11m[q_]:=Module[{f,f1,f2,xt,yt},f1=If[TS==1,f2=SigmaxxYq[q,xt,yt];xt=Re[dd[q]Cosh[Xi]];yt=Im[dd[q]Cosh[Xi]];
f2,0];f=Re[AlphaFar[q]+Alpha0[q]]+Re[f1]]

Sigma12m[q_]:=Module[{f,f1,f2,xt,yt},f1=If[TS==1,f2=SigmaxyYq[q,xt,yt];xt=Re[dd[q]Cosh[Xi]];yt=Im[dd[q]Cosh[Xi]];
f2,0];
f=-Im[AlphaFar[q]+Alpha0[q]]+Re[f1]]

Sigma22m[q_]:=Module[{f,f1,f2,xt,yt},f1=If[TS==1,f2=SigmayyYq[q,xt,yt];xt=Re[dd[q]Cosh[Xi]];yt=Im[dd[q]Cosh[Xi]];
f2,0];
f=BetaFar[q]+Beta0[q]-Re[Alpha0[q]]-Re[AlphaFar[q]]+Re[f1]]

(*Sigma11mlist:=Module[{f},f=Table[Sigma11m[q],{q,ntot}]]
Sigma12mlist:=Module[{f},f=Table[Sigma12m[q],{q,ntot}]]
Sigma22mlist:=Module[{f},f=Table[Sigma22m[q],{q,ntot}]]*)



(*Rotated stress fields in the global z coordinate system*)
Clear[Sigma11z,Sigma22z,Sigma12z]
Clear[zp,zq,Xi,z]
Sigma11z[i_]:=Module[{f,f1,f2,f3},f1=If[TS==1,f2=SxxFinalTemp[i,xt,yt];xt=Re[dd[i]Cosh[Xi]];yt=Im[dd[i]Cosh[Xi]];
f2,0];f3=Sigma11mStresslist[[i]] Cos[Theta[i]]^2-2 Sigma12mStresslist[[i]] Cos[Theta[i]] Sin[Theta[i]]+Sigma22mStresslist[[i]] Sin[Theta[i]]^2;
f=f3+Re[f1]];

Sigma12z[i_]:=Module[{f,f1,f2,f3},f1=If[TS==1,f2=TxyFinalTemp[i,xt,yt];xt=Re[dd[i]Cosh[Xi]];yt=Im[dd[i]Cosh[Xi]];
f2,0];f3=Sigma12mStresslist[[i]] Cos[2 Theta[i]]+(Sigma11mStresslist[[i]]-Sigma22mStresslist[[i]]) Cos[Theta[i]] Sin[Theta[i]];
f=f3+Re[f1]];

Sigma22z[i_]:=Module[{f,f1,f2,f3},f1=If[TS==1,f2=SyyFinalTemp[i,xt,yt];xt=Re[dd[i]Cosh[Xi]];yt=Im[dd[i]Cosh[Xi]];
f2,0];f3=Sigma22mStresslist[[i]] Cos[Theta[i]]^2+Sigma11mStresslist[[i]] Sin[Theta[i]]^2+Sigma12mStresslist[[i]] Sin[2 Theta[i]];
f=f3+Re[f1]];



Clear[zp,zq,Xi,z]
Sigma11zStress[i_]:=Module[{f},f=Sigma11mStresslist[[i]] Cos[Theta[i]]^2-2 Sigma12mStresslist[[i]] Cos[Theta[i]] Sin[Theta[i]]+Sigma22mStresslist[[i]] Sin[Theta[i]]^2];
Sigma12zStress[i_]:=Module[{f},f=Sigma12mStresslist[[i]] Cos[2 Theta[i]]+(Sigma11mStresslist[[i]]-Sigma22mStresslist[[i]]) Cos[Theta[i]] Sin[Theta[i]]];
Sigma22zStress[i_]:=Module[{f},f=Sigma22mStresslist[[i]] Cos[Theta[i]]^2+Sigma11mStresslist[[i]] Sin[Theta[i]]^2+Sigma12mStresslist[[i]] Sin[2 Theta[i]]];


(*Sigma1111z[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi];
f1=Sigma11z[i];
f2=f1;
Xi=XiQ[i];
z=x+I y;
f3=f2]

Sigma1211z[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma12z[i];
Xi=XiQ[i];
z=x+I y;
f2=f1]

Sigma2211z[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma22z[i];
Xi=XiQ[i];
f2=f1;
z=x+I*y;
f3=f2]*)

Sigma1111z[i_]:=Module[{f,f1,f2,f3,f4},Clear[x,y,zp,zq,Xi];
f1=Sigma11zStress[i];
f2=f1;
Xi=XiQ[i];
z=x+I y;
f4=SxxFinalz[x,y];
f3=f2+f4]

Sigma1211z[i_]:=Module[{f,f1,f2,f3,f4},Clear[x,y,zp,zq,Xi,z];
f1=Sigma12zStress[i];
Xi=XiQ[i];
z=x+I y;
f4=TxyFinalz[x,y];
f2=f1+f4]

Sigma2211z[i_]:=Module[{f,f1,f2,f3,f4},Clear[x,y,zp,zq,Xi,z];
f1=Sigma22zStress[i];
Xi=XiQ[i];
f2=f1;
z=x+I*y;
f4=SyyFinalz[x,y];
f3=f2+f4]

Sigma2211ztest[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma22z[i];
Eta=0;
Mu=0;
Xi=XiQ[i];
f2=f1;
z=x+I y;
f3=f2]

PlotSigma11zQ[i_]:=Block[{r,rr},Clear[eta,ss,r1,\[Xi],r];
r=Sigma11z[i];
Xi=Zeta0[i]+I eta;
rr=Plot[r,{eta,0,2Pi},PlotRange->Full]]


Sigma1111zStress[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi];
f1=Sigma11zStress[i];
f2=f1;
Xi=XiQ[i];
z=x+I y;
f3=f2]

Sigma1211zStress[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma12zStress[i];
Xi=XiQ[i];
z=x+I y;
f2=f1]

Sigma2211zStress[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma22zStress[i];
Xi=XiQ[i];
f2=f1;
z=x+I*y;
f3=f2]


PlotSigma12zQ[i_]:=Block[{r,rr},Clear[eta,ss,r1,\[Xi],r];
r=Sigma12z[i];
Xi=Zeta0[i]+I eta;
rr=Plot[r,{eta,0,2Pi},PlotRange->Full]]

PlotSigma22zQ[i_]:=Block[{r,rr},Clear[eta,ss,r1,\[Xi],r];
r=Sigma22z[i];
Xi=Zeta0[i]+I eta;
rr=Plot[r,{eta,0,2Pi},PlotRange->Full]]


PlotSigma11z[i_,xmin_,xmax_,ymin_,ymax_]:=Block[{r0,r},Clear[z,\[Xi]];Unprotect[x,y];
z=x+I y;
r0=Sigma1111z[i];
r=ContourPlot[r0,{x,xmin,xmax},{y,ymin,ymax},FrameLabel->{{"y",""},{"x",""}},ContourLabels->True,PlotLegends->BarLegend[Automatic,All],PlotLabel->"Stress Field \[Sigma]11",Exclusions->None]]

PlotSigma12z[i_,xmin_,xmax_,ymin_,ymax_]:=Block[{r,rr},Clear[r0,u,r0,z,\[Xi]];Unprotect[x,y];
z=x+I y;
r0=Sigma1211z[i];
ContourPlot[r0,{x,xmin,xmax},{y,ymin,ymax},FrameLabel->{{"y",""},{"x",""}},ContourLabels->True,PlotLegends->BarLegend[Automatic,All],PlotLabel->"Stress Field \[Sigma]12",Exclusions->None]]

PlotSigma22z[i_,xmin_,xmax_,ymin_,ymax_]:=Block[{r,rr},Clear[r0,u,r0,z,\[Xi]];Unprotect[x,y];
z=x+I y;
r0=Sigma2211z[i];
ContourPlot[r0,{x,xmin,xmax},{y,ymin,ymax},FrameLabel->{{"y",""},{"x",""}},ContourLabels->True,PlotLegends->BarLegend[Automatic,All],PlotLabel->"Stress Field \[Sigma]22",Exclusions->None]]


End[] 
Protect@@Names["StressFields`*"]
EndPackage[] 

BeginPackage["StressFields`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","CombinDCDTEMP`"}] 
EndPackage[] 



















