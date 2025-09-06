(* ::Package:: *)

BeginPackage["run0`"]
Unprotect@@Names["run0`*"];
ClearAll@@Names["run0`*"];

(*In this file stress fields of the system are calculated using the multipole 
method. Then, by using the correct coefficients and removing regular fields 
(by equating Eta_nmpq and Mu_nmpq ), stress fields can be calculated at any point 
in the model.*)

SigmaEtaEta0::usage=" \[Sigma]\[Eta]\[Eta]; stress field in yq local coordinate system"; 
Sigma11m0::usage=" sigma11; stress field in yq local coordinate system"; 
Sigma12m0::usage=" sigma12; stress field in yq local coordinate system"; 
Sigma22m0::usage=" sigma22; stress field in yq local coordinate system"; 
Sigma11z0::usage=" sigma11; stress field in zq local coordinate system"; 
Sigma12z0::usage=" sigma12; stress field in zq local coordinate system"; 
Sigma22z0::usage=" sigma22; stress field in zq local coordinate system"; 
Sigma1111z0::usage=" sigma11; stress field in z global coordinate system"; 
Sigma1211z0::usage=" sigma12; stress field in z global coordinate system"; 
Sigma2211z0::usage=" sigma22; stress field in z global coordinate system"; 

Begin["`Private`"] 



AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];

Needs["SubKernels`LocalKernels`"];
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]]);




Get["run1.m"];
Needs["run1`"];
CloseKernels[subkernel]; initsub[];
(*Get["SIF.m"];*)
Get["linearSolve.m"];
Get["DisplacementFields.m"];
Needs["DisplacementFields`"];

Print["Finite Boundary Calculation. run0 start"];




Eta1[n_,m_,p_,q_]:=Module[{f1},f1=If[Abs[zpq[p,q]]>=Abs[ltest],0,0]];

Eta2:=Module[{f},f=Table[Eta1[i,m,p,q],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}]];
Eta3=Eta2;

Mu1[n_,m_,p_,q_]:=Module[{f1},f1=If[Abs[zpq[p,q]]>=Abs[ltest],0,0]];

Mu2:=Module[{f},f=Table[Mu1[i,m,p,q],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}]];
Mu3=Mu2;

Eta3>>Eta.2;
Mu3>>Mu.2;



Unprotect[Eta,Mu];
Eta = <<Eta.2;
Mu = <<Mu.2;





(*Stresses Caused by the qth inclusion in the local yq coordinate system.*)
SigmaEtaEta0[q_]:=Block[{f},f=0.5*(Beta0[q]+BetaFar[q]-Re[Tau0[q]+TauFar[q]])];
Unprotect[z,Xi,x,y];
(*z = x + I y*)
Sigma11m0[q_]:=Module[{f},f=Re[AlphaFar[q]+Alpha0[q]]];
Sigma12m0[q_]:=Module[{f},f=-Im[AlphaFar[q]+Alpha0[q]]];

Sigma11mlist:=Module[{f},f=Table[Sigma11m0[q],{q,ntot}]];
Sigma12mlist:=Module[{f},f=Table[Sigma12m0[q],{q,ntot}]];

Sigma22m0[q_]:=Module[{f},f=BetaFar[q]+Beta0[q]-Re[Alpha0[q]]-Re[AlphaFar[q]]];
Sigma22mlist:=Module[{f},f=Table[Sigma22m0[q],{q,ntot}]];

(*Rotated stress fields in the global z coordinate system*)
Clear[Sigma11z0,Sigma22z0,Sigma12z0];


(*Stresses Caused by the qth inclusion in the local zq coordinate system.*)
Clear[zp,zq,Xi,z];
Sigma11z0[i_]:=Module[{f},f=Sigma11mlist[[i]] Cos[Theta[i]]^2-2 Sigma12mlist[[i]] Cos[Theta[i]] Sin[Theta[i]]+Sigma22mlist[[i]] Sin[Theta[i]]^2];
Sigma12z0[i_]:=Module[{f},f=Sigma12mlist[[i]] Cos[2 Theta[i]]+(Sigma11mlist[[i]]-Sigma22mlist[[i]]) Cos[Theta[i]] Sin[Theta[i]]];
Sigma22z0[i_]:=Module[{f},f=Sigma22mlist[[i]] Cos[Theta[i]]^2+Sigma11mlist[[i]] Sin[Theta[i]]^2+Sigma12mlist[[i]] Sin[2 Theta[i]]];


(*Stresses Caused by the qth inclusion in the global z coordinate system.*)
Sigma1111z0[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi];
f1=Sigma11z0[i];
f2=f1;
Xi=XiQ[i];
z=x+I y;
f3=Re[f2]];

Sigma1211z0[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma12z0[i];
Xi=XiQ[i];
z=x+I y;
f2=Re[f1]];

Sigma2211z0[i_]:=Module[{f,f1,f2,f3},Clear[x,y,zp,zq,Xi,z];
f1=Sigma22z0[i];
Xi=XiQ[i];
f2=f1;
z=x+I y;
f3=f2];

Print["Finite Boundary Calculation. run0 END"];


End[] 
Protect@@Names["run0`*"]
EndPackage[] 

(*BeginPackage["run0`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","PackageManipulations`","run1`","CombinMMBE`","eta0p`","AddInclusion`","Addone`"}] 
EndPackage[] *)
BeginPackage["run0`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`"}] 
EndPackage[] 


































