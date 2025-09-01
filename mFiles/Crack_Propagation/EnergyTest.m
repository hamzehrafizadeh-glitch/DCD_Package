(* ::Package:: *)

BeginPackage["EnergyTest`"]
Unprotect@@Names["EnergyTest`*"];
ClearAll@@Names["EnergyTest`*"];


(*here I imagine that in the system under calculation,a crack path consists of two micro-cracks.*)

(*At this point and before determining whether the crack can propagate,deflect,or arrest,
it should be determined where is the crack tip and what is the corresponding fracture energy.*)



g::usage="Maximum Energy release rate of a propagating crack.:";
PropagationTest::usage="Check if the crack energetically can propagate?:";
Energy::usage="";
gt1::usage="";
AppendCrackEnergy::usage="AppendCrackEnergy[p_]";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
(*Get["eta0p.m"]*)
Needs["eta0p`"];
Needs["AddInclusion`"]

Unprotect[ninitial,lltest]
ninitial=readstart[[1]]
lltest=readstart[[2]]
(*I have to find a way for avoiding or adding ttest in the code. ttest is used in the run4 file.*)


g[p_, ttest_] := 
  Block[{fg, kk, k2, f,ii}, 
ii=2;
(*Print["ntot for energy =   ",  p];*)
fg = (Chi0 + 1)/8/Mu0;(*1/GPa*)
   kk = If[ntot==0||ntot-1 == ninitial, SIF`kcrack[p, p], 
     SIF`kcrack[p - ttest+ii, p - ttest+ii-1 ]
    (*SIF`kcrack[p - ttest, p - ttest + 1]*)];
   k2 = kk*Conjugate[kk];
   f = Re[1000*fg*k2]];
(*At this point and before determining whether the crack can propagate,deflect,or arrest,it should be determined where is the crack tip and what is the corresponding fracture energy.*)


(*l1temp1[[i]] is defined in cracklength function in run4*)
(*j will be used regarding run4 file.*) 
PropagationTest[p_,j_, ttest_] := 
 Block[{f,f1,kk,g1,g2,i,ninitial1,ii}, 
ii=2;
Print["ntot  ", ntot];
ninitial1=If[ ninitial==0,1,ninitial];
f1 = If[j==1,EndPoint[p],EndPoint[p-ttest+ii]];(*Sum[l1temp[[i]], {i, ninitial1, Times @@ Dimensions[l1temp]}];*)
g1=If[j==1,G0,gc];g2=g[p, ttest];
Print["gc = ", g1]; Print["g[ntot] = ", g2];
kk = If[ntot==1||ntot-1 == ninitial, SIF`kcrack[p, p], 
     SIF`kcrack[p - ttest+ii, p - ttest+ii-1]];
Print["SIF =  ", kk];
PutAppend[List[Re[kk],Im[kk],g1,g2,f1,e[p-ttest+ii]],FileNameJoin[{FilePath, "output", "crackenergy"}]];
  f = If[g2 < g1, Print["Crack is arrested!"];, 
    Print["Crack is propagating!"]]];




(*Check these functions*)
Energy[p_, q_] := Block[{x, y, r, coeff, field, f, eta},
  x = Dq[p] Cos[eta] Cosh[zeta];
  y = Dq[p] Sin[eta] Sinh[zeta];
  e1 = SigmaEtaEta[p]^2;
  zeta = Zeta0[p]; eta = eta0[p];
  r = Sqrt[x^2.0 + y^2.0];
  Print["radius  ", r];
  coeff = 1000(Pi*r/4.0/Mu0)*(Chi0 + 1);(*1/GPa*)
  
  Print["Coeff ", coeff];
  ave = If[p == q, 0, SEtaEtaAve[p, q, zeta, eta]^2.0];
  Xi = zeta + I*eta;
  field = e1 + ave;
  f = coeff*field]

gt[p_, ttest_] := Block[{fg, kk, k2, f}, Print["ntot =   ", p];
  kk = If[ntot - 1 == ninitial, Energy[p, p], 
    Energy[p - ttest + 1, p - ttest]]]


(*Write information into the crackenergy file.*)
AppendCrackEnergy[j_,p_,ttest_]:=Block[{f1,g1,g2,kk,fg,k2,f,i},
f1 = EndPoint[p];Sum[l1temp[[i]], {i, ninitial, Times @@ Dimensions[l1temp]}];
g1=GQ[j];
Print["gc = ", g1]; 
kk = StoredStressInc[p,ttest];
Print["SIF =   ", kk];
fg = (Chi0 + 1)/8/Mu0;(*1/GPa*)
k2 = kk*Conjugate[kk];
  g2 = Re[1000*fg*k2];Print["g[ntot] = ", g2];
PutAppend[List[Re[kk],Im[kk],g1,g2,f1,e[p],p,ntot,AppendCrack],FileNameJoin[{FilePath, "output", "crackenergy"}]];
  f = If[g2 < g1, Print["Crack is arrested!"];, 
    Print["Crack is propagating!"]]];


End[] 
Protect@@Names["EnergyTest`*"]
EndPackage[] 

BeginPackage["EnergyTest`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","AddInclusion`","Addone`","InterfaceFunc`"}] 
EndPackage[] 








































