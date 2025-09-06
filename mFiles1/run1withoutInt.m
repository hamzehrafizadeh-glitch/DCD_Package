(* ::Package:: *)

BeginPackage["run1`"]
Unprotect@@Names["run1`*"];
ClearAll@@Names["run1`*"];

Begin["`Private`"] 



(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"]





Eta1[n_,m_,p_,q_]:=Module[{f1},f1=If[Abs[zpq[p,q]]>=Abs[ltest],EtaFarPQ[n,m,p,q],EtaClosePQ[n,m,p,q]]]

Eta2:=Module[{f},f=Table[Eta1[i,m,p,q],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}]]
Eta3=Eta2;
Eta4=Table[If[Abs[Eta3[[i,m,p,q]]]<=10^-4,0,Eta3[[i,m,p,q]]],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}];

Mu1[n_,m_,p_,q_]:=Module[{f1},f1=If[Abs[zpq[p,q]]>=Abs[ltest],MuFarPQ[n,m,p,q],MuClosePQ[n,m,p,q]]]

Mu2:=Module[{f},f=Table[Mu1[i,m,p,q],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}]]
Mu3=Mu2;
Mu4=Table[If[Abs[Mu3[[i,m,p,q]]]<=10^-4,0,Mu3[[i,m,p,q]]],{i,n+1},{m,nmax1},{p,ntot},{q,ntot}];

Eta4>>Eta.1;
Mu4>>Mu.1;

(*Quit[]*)







End[] 
Protect@@Names["run1`*"]
EndPackage[] 

BeginPackage["run1`",{"geometry`", "ExpansionCoefficient`" }] 
EndPackage[] 






















