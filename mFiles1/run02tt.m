(* ::Package:: *)

BeginPackage["run02`"]
Unprotect@@Names["run02`*"];
ClearAll@@Names["run02`*"];


(*This file calculates the stress fields caused by the inclusions on the points 
of the mesh of the grain boundary. This information (SMM) will be used as the 
supplementary boundary condition for BE method.*)


(*The static analysis can be finalised in six steps. 
Calculating the influence coefficient matrices (step 3) and the stress and displacement 
fields (step 6) are the most computationally expensive steps in the analysis. In the 
crack propagation problems, where an iterative solution process is required, the 
iteration is performed over steps 4 and 5.
*)


Tnfield::usage ="Tnfield[t]: Tangential stress field for element with t Degree rotation";

Tsfield::usage ="Tsfield[t]: Normal stress field for element with t Degree rotation";


Begin["`Private`"] 





AppendTo[$Path, 
  ToFileName[{$HomeDirectory, 
    "Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];
(*Get["PackageManipulations.m"]
Needs["PackageManipulations`"]*)
Get["GenerateMesh.m"];
Needs["GenerateMesh`"];



Needs["SubKernels`LocalKernels`"];
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]]);
Print["run02 Calculation."];




Get["run0.m"];
Needs["run0`"];


(*points with fixed normal displacements
Points have to be written as a string like {2,4,7}*)

Clear[Sn]
Sn[j_]:=Block[{f,jj,j1},jj=j;
test = Table[1, {i, NTOTAL}];
Do[test[[j1]]=0,{j1,jj}];
test
]

(*points with fixed tangential displacements*)

Clear[Ss]
Ss[j_]:=Block[{f,jj,j1},jj=j;
test = Table[1, {i, NTOTAL}];
Do[test[[j1]]=0,{j1,jj}];
test
]

(*Choose j*)

(*jn={NTOTAL-1,NTOTAL}
js={NTOTAL-1}*)

jn={2n1+n2+1,2n1+2n2}
js={2n1+2n2}
(*js={};{2n1+2n2}*)
(**)


Snn=Sn[jn];
Sss=Ss[js];


Clear[Tn,Ts]
Tn[t_,i_]:=Block[{f,S11,S22,S12},f=Simplify[(S11+S22)/2. +Cos[2*t](S11-S22)/2.  +S12*Sin[2*t]];
S11:=Sigma1111z0t[i]-s11;
S22:=Sigma2211z0t[i]-s22;
S12:=Sigma1211z0t[i]-s12;
f]

Ts[t_,i_]:=Block[{f,S11,S22,S12},f=Simplify[-Sin[2*t](S11-S22)/2.  +S12*Cos[2*t]];
S11:=Sigma1111z0t[i]-s11;
S22:=Sigma2211z0t[i]-s22;
S12:=Sigma1211z0t[i]-s12;
f]

Clear[Tnfield,Tsfield]
Tnfield[t_]:=Block[{f,S11,S22,S12},f=Simplify[(S11+S22)/2. +Cos[2*t](S11-S22)/2.  +S12*Sin[2*t]];
S11:=s11;
S22:=s22;
S12:=s12;
f]

Tsfield[t_]:=Block[{f,S11,S22,S12},f=Simplify[-Sin[2*t](S11-S22)/2.  +S12*Cos[2*t]];
S11:=s11;
S22:=s22;
S12:=s12;
f]


SnMM := Block[{r, f,i,j,t}, r :=Simplify[ Sum[Tn[t,i], {i, 1, ntot}]+Tnfield[t]];
  f = Table[x = Re[point[[j]]]; y = Im[point[[j]]];t=-XY[[j,4]]+Pi/2; If[Snn[[j]]==1, r,0], {j, 1, NTOTAL}]]

SsMM := Block[{r, f,i,j,t}, r :=Simplify[ Sum[Ts[t,i], {i, 1, ntot}]+Tsfield[t]];
  f = Table[x = Re[point[[j]]]; y = Im[point[[j]]];t=-XY[[j,4]]+Pi/2; If[Sss[[j]]==1, r,0], {j, 1, NTOTAL}]]

Print["XY"];
XY = << XY.dat;
(*Print["XY", XY];*)
point = Block[{i},Table[XY[[i, 1]] + I* XY[[i, 2]], {i, 1, NTOTAL}]];

SnMM;
SsMM ;

(*SnMM=SnMM1 Snn;
SsMM=SsMM1 Sss ;*)



Print["Matrixes SMM Start."];
Print["Time"];
Print[DateList[]];
Table[Sigma1111z0t[i]=Sigma1111z0[i],{i,1,ntot}];
Table[Sigma2211z0t[i]=Sigma2211z0[i],{i,1,ntot}];
Table[Sigma1211z0t[i]=Sigma1211z0[i],{i,1,ntot}];

Print["Matrixes SMM End."];
Print["Time"];
Print[DateList[]];

SnMM;
SsMM ;

(*SnMM=SnMM1 Snn;
SsMM=SsMM1 Sss ;*)

SSMM=Chop[Join[SnMM, SsMM]];



(*Displacement fixed for points j\[Equal]2n1+n2+1||j\[Equal]2n1+2n2*)
(*If pure tension, then SUMM=0*)

Clear[Un,Us]
Clear[uxzt,uyzt]
Un[t_,i_]:=Block[{f,uxt,uyt},f=Simplify[uxt Cos[t]+uyt Sin[t]];
uxt:=uxzt[i];
uyt:=uyzt[i];
f]

Us[t_,i_]:=Block[{f,uxt,uyt},f=Simplify[uyt Cos[t]-uxt Sin[t]];
uxt:=uxzt[i];
uyt:=uyzt[i];
f]



Clear[SUxMM,SUyMM,SUMM]
SUnMM := Block[{r, f,i,j,t}, r :=Simplify[ Sum[Un[t,i], {i, 1, ntot}]];
  f = Table[x = Re[point[[j]]]; y = Im[point[[j]]];t=-XY[[j,4]]+Pi/2; If[Snn[[j]]==1, 0,r], {j, 1, NTOTAL}]]

SUsMM := Block[{r, f,i,j,t}, r :=Simplify[ Sum[Us[t,i], {i, 1, ntot}]];
  f = Table[x = Re[point[[j]]]; y = Im[point[[j]]];t=-XY[[j,4]]+Pi/2; If[Sss[[j]]==1, 0,r], {j, 1, NTOTAL}]]


SUnMM ;
SUsMM;


Table[uxzt[i]=uxz[i],{i,1,ntot}];
Table[uyzt[i]=uyz[i],{i,1,ntot}];

SUnMM;
SUsMM;

(*SUnMM=SUnMM1 (1-Snn);
SUsMM=SUsMM1 (1-Sss);*)

SUMM = Chop[Join[SUnMM, SUsMM]];



SMM = SUMM +SSMM; 
Print["Stress fields at the boundary"]
(*In the version of mathematica it needs Chop!*)
(*SMM1r=Table[If[Abs[Re[SMM[[i]]]]\[LessEqual]10^-9,0,Re[SMM[[i]]]],{i, 1, 2*NTOTAL}];
SMM1im=Table[If[Abs[Im[SMM[[i]]]]\[LessEqual]10^-9,0,Im[SMM[[i]]]],{i, 1, 2*NTOTAL}];
SMM1=Table[SMM1r[[i]]+I*SMM1im[[i]],{i, 1, 2*NTOTAL}]*)
SMM >> SMM.dat;
Print["SMM",SMM];

Print["SUMM",SUMM];
Print["SSMM",SSMM];

ClearSystemCache[];
 CloseKernels[subkernel]; initsub[];



(*Solve the BE code*)



CloseKernels[subkernel]; initsub[]; 
Print["Fieldsxy Statrt"];
Get["Fieldsxy.m"];
Needs["Fieldsxy`"];
Print["Fieldsxy End"];
(*PackageReload["Fieldsxy`", KillShadowing -> True]; *)
CloseKernels[subkernel]; initsub[]; 



(*Change expansion coefficient to the original values and combine BE and MM by 
loading CombinMMBE file.*)


(*Unprotect[Eta, Mu]*)
ClearAll[Eta, Mu];

Eta = << Eta.1;
Mu = << Mu.1;




Get["CombinMMBE.m"];
Needs["CombinMMBE`"];


Print["run02 Calculation.END"];


End[] 
Protect@@Names["run02`*"]
EndPackage[] 

BeginPackage["run02`",{"geometry`", "run1`","run0`","Fieldsxy`","CombinMMBE`" (*,"eta0p`","AddInclusion`","Addone`"*)}] 
EndPackage[] 


(*BeginPackage["run02`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","PackageManipulations`","run1`","run0`","Fieldsxy`","CombinMMBE`" (*,"eta0p`","AddInclusion`","Addone`"*)}] 
EndPackage[] *)












































































