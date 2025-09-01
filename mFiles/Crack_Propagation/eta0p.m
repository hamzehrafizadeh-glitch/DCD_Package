(* ::Package:: *)

BeginPackage["eta0p`"]
Unprotect@@Names["eta0p`*"];
ClearAll@@Names["eta0p`*"];




eta0::usage="eta0[p]: It determines an angle, in which the tangential stress on the surface of the inclusion q is maximum. "; 
myArcTan::usage="myArcTan[x_,y_], it returns the angle in a range of [0,2Pi]. ";
readstart::usage = "It is a function that opens and reads the start file.";
writestart::usage = "It is a function that opens and writes the start file.";
ninitial::usage = "Start file first row.";
lltest::usage = "start file second row.";
eta0L::usage="";
SIFTemp::usage="SIFTemp[p]: It gives SIF associated to the Temperature variation at the crack tip of crack p";




Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)

(*Get["PackageManipulations.m"];
Needs["PackageManipulations`"];*)

Get["geometry.m"]
Needs["geometry`"]
If[BE==1,
Print["BE, I am in eta0 file"];
(*Get["run02.m"];*)
Needs["run02`"];
(*Get["run02.m"];
Needs["run02`"];*)
,
Print[" I am in eta0 file"];
(*Get["run1.m"];*)
Needs["run1`"];
CloseKernels[subkernel]; initsub[];
Needs["linearSolve`"]];


Get["DisplacementFields.m"];
Needs["DisplacementFields`"];
Needs["CombinDCDTEMP`"];

(*PackageReload["run1`",KillShadowing\[Rule]True];
CloseKernels[subkernel];initsub[];
PackageReload["SIF`",KillShadowing\[Rule]True];
CloseKernels[subkernel];initsub[];
PackageReload["eta0p`",KillShadowing\[Rule]True];*)


(*myArcTan[x_,y_] := Piecewise[{{Pi + ArcTan[y/x], x < 0}}, ArcTan[y/x]];*)
myArcTan[x_,y_] := Piecewise[{{2.0*Pi + ArcTan[x,y], (x < 0.0 && y<0)}}, ArcTan[x,y]];


writestart:=Block[{s8,rtest},s8 = OpenWrite[StartFile];
Write[s8,ninitial];
Write[s8,lltest];
Close[s8]]

readstart:=Block[{s8,rtest},s8 = OpenRead[StartFile];
rtest =Flatten[ Table[ReadList[s8, Expression, 1], {j, 1+1}]];
Close[s8];
rtest]

ninitial=readstart[[1]];
lltest=readstart[[2]];

(*Note: inc<=ttest-2*)(*decide how many inclusions left at the end of the crack path to 
predict path more accurate.*)





eta0test[p_]:=Block[{r, r1, r2,r3,t,f1,f,eta1},Unprotect[eta,eta1];Clear[Xi];
If[e[p]<=0.01,t=100,t=2];
 r1 = SigmaEtaEta[p];
 Xi = 1.*Zeta0[p] + I*eta;
 r2 = SEtaEta[p]; 
eta1=eta;
 r3 = r1 + r2/.{Conjugate[eta] -> eta};
 r=If[BE==1,Re[r3],Re[r1]];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2, 
    Last[NMaximize[{r, -Pi/t < eta < Pi/t}, eta, AccuracyGoal -> 5]], 
     Last[NMaximize[{r, Pi-Pi/t < eta < Pi+Pi/t}, eta, AccuracyGoal -> 5]]];
f=eta /.f1]



eta0t[p_]:=Block[{r, r1, r2,r3,t,f1,f,eta1,eta2},Unprotect[eta,eta1,eta2];Clear[Xi];
If[e[p]<=0.01,t=100,t=20];
 r1SF = SigmaEtaEta[p];
 Xi = 1.5*Zeta0[p] + I*eta;
r2 =If[BE==1, SEtaEta[p],0]; 
eta1=eta;
r2TS =If[TS==1, SEtaEtaTS[p],0]; 
eta2=eta;
r1=r1SF+r2TS;
 r3 = (r1 + r2)/.{Conjugate[eta] -> eta};
 r=Re[r3];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2, 
    Last[NMaximize[{r, -Pi/t < eta < Pi/t}, eta, AccuracyGoal -> 5]], 
     Last[NMaximize[{r, Pi-Pi/t < eta < Pi+Pi/t}, eta, AccuracyGoal -> 5]]];
f=eta /.f1]


eta0[p_]:=Block[{r, r1, r2,r3,t,f1,f,eta1,eta2,t1,r1SF,r2TS},Unprotect[eta,eta1,eta2];Clear[Xi];t1=2.;Unprotect[dd];
(*Clear[dd];*)dd[p] = t1*Sqrt[l1[p]^2 - l2[p]^2] Exp[I Theta[p]];
If[e[p]<=0.01,t=100,t=20];
 r1SF = SigmaEtaEta[p];
 Xi = 1.3*Zeta0[p] + I*eta;
r2 =If[BE==1, Print["eta0, BE"];SEtaEta[p],0]; 
eta1=eta;
r2TS =If[TS==1, SEtaEtaTS[p],0]; 
eta2=eta;
r1=r1SF+r2TS;
 r3 = (r1 + r2)/.{Conjugate[eta] -> eta};
 r=Re[r3];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2, 
    Last[NMaximize[{r, -Pi/t < eta < Pi/t}, eta, AccuracyGoal -> 5]], 
     Last[NMaximize[{r, Pi-Pi/t < eta < Pi+Pi/t}, eta, AccuracyGoal -> 5]]];
(*Clear[dd];*)Get["geometry.m"];
Needs["geometry`"];
f=eta /.f1]


eta0L[p_]:=Block[{r, r1, r2,r3,t,f1,f,eta1},Unprotect[eta];Clear[Xi];
If[e[p]<=.01,t=100,t=2];
 r1 = SigmaEtaEta[p];
 Xi = 1.*Zeta0[p] + I*eta;
r2 = SEtaEta[p]; 
eta1=eta;
 r3 = r1 + r2/.{Conjugate[eta] -> eta};
 r=If[BE==1,Re[r3],Re[r1]];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2,  
     Last[NMaximize[{r, Pi-Pi/t < eta < Pi+Pi/t}, eta, AccuracyGoal -> 5]],
     Last[NMaximize[{r, -Pi/t < eta < Pi/t}, eta, AccuracyGoal -> 5]]];
f=eta /.f1]


(*eta0[p_]:=Block[{r, r1, r2,r3,rt1,rt2,rt3,rrt,f1,f},Unprotect[eta];Clear[Xi];
 r1 = SigmaEtaEta[p];
 rt1=Tau0[p]+TauFar[p];
 Xi = 1.*Zeta0[p] + I*eta;
 rt2=Im[rt1];
 rrt=r1-rt2;
 rt3=Re[rrt*Conjugate[rrt]]/.{Conjugate[eta] -> eta};

 r2 = SEtaEta[p]; 
 r3 = r1 + r2/.{Conjugate[eta] -> eta};
 r=If[BE==1,Re[r3],Re[rt3]];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2, 
    Last[NMaximize[{r, -Pi/4 < eta < Pi/4}, eta, AccuracyGoal -> 5]], 
     Last[NMaximize[{r, 3Pi/4 < eta < 5*Pi/4}, eta, AccuracyGoal -> 5]]];
f=eta /.f1]*)


eta0SIFt[p_]:=Block[{f,r1,r2,k11,k22},If[ntot==1,q=p,q=p-1];
k11=Re[kcrack[p,q]];
k22=Im[kcrack[p,q]];
r1=ArcCos[(3*k22^2-Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
r2=ArcCos[(3*k22^2+Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
f=If[k22>0,-r2,r2]]
(*For this one I have to write a new NewCenter[p_] one considering our new angle*)





SIFTemp1[p_]:= Block[{etatemp,r,r12,f22,f12,f, Xi, z1,x1,y1,ztemp,slop,sloptemp},
  etatemp = eta0[p];
  ztemp = Dq[p] Cosh[Zeta0[p] + I*etatemp];
 r = Re[SigmayyYq[p,x1,y1]];
  r12=Re[SigmaxyYq[p,x1,y1]];
  y1 = Im[ztemp] ;
  x1=Re[ztemp];
  f22= 
   Sqrt[(-Dq[p]+ Re[ztemp])*2.Pi]r;
  f12=
   Sqrt[(-Dq[p]+ Re[ztemp])*2.Pi]r12;
	f=f22+I*f12
  ]

SIFTemp1[i_]:=Block[{f,xini,xtemp,ytemp},xini=Dq[i];f=Simplify[Sqrt[2Pi*((xtemp-Re[xini])^2 +(ytemp-Im[xini])^2)](SigmayyYq[i,xtemp,ytemp]+I*SigmaxyYq[i,xtemp,ytemp])];
etatemp = eta0[i];
  ztemp = Dq[i] Cosh[Zeta0[i] + I*etatemp];
xtemp=Re[ztemp];
ytemp=Im[ztemp];
f]

SIFTemp[i_]:=Block[{f,xini,xtemp,ytemp},xini=Dq[i];f=Simplify[Sqrt[2Pi*((xtemp-Re[xini])^2 +(ytemp-Im[xini])^2)](SigmayyYq[i,xtemp,ytemp]+I*SigmaxyYq[i,xtemp,ytemp])];
xtemp=Re[(EndPoint[i]-zcenter[i])Exp[-I*Thetaprim[i]]];
ytemp=Im[(EndPoint[i]-zcenter[i])Exp[-I*Thetaprim[i]]];
f]


(*By selecting r2 as a propagation angle,the crack propagation is confined to merely 
the RHS propagation.This means that the angle is between [-Pi/2 ,Pi/2]*)

(*In order to insure that the energy release associated with the crack extension is 
maximum,the sign of the propagation angle $\theta$ should be opposite to the sign of
$K_{II}$ \cite{wei1982nonlinear}*)




End[] 
Protect@@Names["eta0p`*"]
EndPackage[] 

BeginPackage["eta0p`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","CombinMMBE`","CombinDCDTEMP`"}] 
EndPackage[] 





















































