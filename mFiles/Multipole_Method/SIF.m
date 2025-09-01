(* ::Package:: *)

      


BeginPackage["SIF`"]
Unprotect@@Names["SIF`*"];
ClearAll@@Names["SIF`*"];


(*In the process of correction of the stress intensity factor for crack 
propagation, I do not need to consider the effect of all inclusions. 
In another word, the stored energy in the inclusion p is a summation 
of the effect of all inclusions. Crack propagation means that the stored 
energy of the inclusion p caused by the inclusion q has to be released, 
while inclusion q is part of the main crack.*)

k::usage=" k[q] = kI[q]+i KII[q];  k shows SIF i and ii (mode I and II) for RHS of the inclusion q "; 
kminus::usage=" kminus[q] = kI[q]+i KII[q];  kminus shows SIF i and ii (mode I and II) for LHS of the inclusion q"; 

means11::usage="means11[p,q]: Is an average stored stress s11 in an inclusion p because of the presence of other inclusion (q) in the system."; 
means12::usage="means12[p,q]: Is an average stored stress s12 in an inclusion p because of the presence of other inclusion (q) in the system."; 
means22::usage="means22[p,q]: Is an average stored stress s22 in an inclusion p because of the presence of other inclusion (q) in the system."; 
SEtaEtaAveG::usage="SEtaEtaAve[p_, q_, zeta_, eta_]: Is an average stored stress SigmaEtaEta in an inclusion p because of the presence of other inclusion (q) in the global coordinatesystem."; 
SEtaEtaAve::usage="SEtaEtaAve[p_, q_, zeta_, eta_]: Is an average stored stress SigmaEtaEta in an inclusion p because of the presence of other inclusion (q) in the pth local coordinate system."; 
kcrack::usage="kcrack[p,q]: is an SIFs of the inclusion p, only if inclusion p is part of the crack propagation path. "; 
kcrackminus::usage="kcrackminus[p,q]: is a LHS SIFs of the inclusion p, only if inclusion p is part of the crack propagation path. "; 
(*StoredStressInc::usage="StoredStressInc[j_,ttest_]";*)
Reload::usage="Reload Package In SIF.m";
eta0SIF::usage="eta0SIF[p]: It determines an angle, in which the tangential stress on the surface of the inclusion q is maximum. STILL IT IS NOT WORKING BE METHOD ";



Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["PackageManipulations.m"];
Needs["PackageManipulations`"];

Print["I am in SIF File."]
Get["DisplacementFields.m"];
Needs["DisplacementFields`"];

Needs["eta0p`"];
Needs["run0`"];



(*Reload := {Get["geometry.m"];PackageReload["run1`", KillShadowing -> True];
   CloseKernels[subkernel]; initsub[];
Print["Rload!!!!!!!"];
   PackageReload["linearSolve`", KillShadowing -> True];
   Needs["SIF`"];
   CloseKernels[subkernel]; initsub[]};*)

Reload := {Get["geometry.m"];If[BE==1,
Print["BE, I am in Reloading files"];
PackageReload["run02`", KillShadowing -> True],PackageReload["run1`", KillShadowing -> True]];
   CloseKernels[subkernel]; initsub[];
Print["Rload!!!!!!!"];
   PackageReload["linearSolve`", KillShadowing -> True];
   Needs["SIF`"];
   CloseKernels[subkernel]; initsub[]};




Clear[k]
SIFTempFile = FileNameJoin[{FilePath, "output", "SIFTemp"}];
k0[q_]:=Module[{f,f1,f2},f1=If[TS==1,SIFTemp[q],0];PutAppend[List[f1],SIFTempFile];
Print["SIFTemp = ",f1];f2=(-2 Mu0 Sqrt[Pi/dd[q]/Exp[-I Theta[q]]] Conjugate[Sum[i (abmnpq`a[i,q]+Conjugate[abmnpq`b[i,q]]),{i,n}]])+f1]
kminus[q_]:=Module[{f,f1,f2},f2=(-2 Mu0 Sqrt[Pi/dd[q]/Exp[-I Theta[q]]] Conjugate[Sum[Power[-1,i-1] i (abmnpq`a[i,q]+Conjugate[abmnpq`b[i,q]]),{i,n}]])]


k1[q_]:=InterfaceFunc`SIFLinemethodLengtht[q,0.5*l1[q]e[q]^2];(**)
(*k1[q_]:=InterfaceFunc`SIFLinemethodLengtht[q,0.7*l1[q]e[q]^2];(*I HAVE CHANGED THIS SIZE*)*)

k[q_]:= If[num2==1,k0[q],k1[q],k1[q]]



(*

(*This fined the effect of particle q on the p one!*)
Clear[mean1,mean2, means11, means22, means12]
mean1[p_, q_] := 
  Module[{f}, 
   f = 16*Mu0*( Re[abmnpq`anpq[1, q, p]]/dd[p])];
mean2[p_, q_] := 
  Module[{f}, 
   f = 8 Mu0*(
        (abmnpq`bnpqminus[1, q, p] - 
          abmnpq`anpq[1, q, p]/v0[q]/v0[q])/dd[p])];

means22[p_, q_] := Module[{f}, f = Re[mean1[p, q] + mean2[p, q]]/2]
means11[p_, q_] := Module[{f}, f = Re[mean1[p, q] - mean2[p, q]]/2]
means12[p_, q_] := Module[{f}, f = Im[mean1[p, q] + mean2[p, q]]/2]

(*This calculates the new SIFs of the inclusion p by taking into \
account the effect of the inclusion q *)

kcrack[p_, q_] := 
 Block[{r1, r2, r3}, r1 = means22[p, q]*Sqrt[Pi*l1[p]];
  r2 = means12[p, q]*Sqrt[Pi*l1[p]];
  r3 = r1 + I*r2 + k[p]]

kcrackminus[p_, q_] := 
 Block[{r1, r2, r3}, r1 = means22[p, q]*Sqrt[Pi*l1[p]];
  r2 = means12[p, q]*Sqrt[Pi*l1[p]];
  r3 = r1 + I*r2 + kminus[p]]*)

Clear[mean1,mean2, means11, means22, means12]
mean1[p_, q_] := 
  Module[{f}, 
   f = 16*Mu0*( Re[abmnpq`anpq[1, q, p]]/dd[p]/Exp[-I Theta[p]])];
mean2[p_, q_] := 
  Module[{f}, 
   f = 8 Mu0*((abmnpq`bnpqminus[1, q, p] - 
          abmnpq`anpq[1, q, p]/v0[q]/v0[q])/dd[p]/Exp[-I Theta[p]])];

means22[p_, q_] := Module[{f}, f = Re[mean1[p, q] + mean2[p, q]]/2]
means11[p_, q_] := Module[{f}, f = Re[mean1[p, q] - mean2[p, q]]/2]
means12[p_, q_] := Module[{f}, f = Im[mean1[p, q] + mean2[p, q]]/2]

(*This calculates the new SIFs of the inclusion p by taking into \
account the effect of the inclusion q *)

kcrack[p_, q_] := 
 Block[{r1, r2, r3}, r1 = means22[p, q]*Sqrt[Pi*l1[p]];
  r2 = means12[p, q]*Sqrt[Pi*l1[p]];
  r3 = r1 + I*r2 + k[p]]

kcrackminus[p_, q_] := 
 Block[{r1, r2, r3}, r1 = means22[p, q]*Sqrt[Pi*l1[p]];
  r2 = means12[p, q]*Sqrt[Pi*l1[p]];
  r3 = r1 + I*r2 + kminus[p]]


(*SEtaEtaAve[p_, q_] := Block[{f},
  x = Re[zcenter[p]] + Dq[p] Cos[eta] Cos[Theta[p]] Cosh[Zeta0[p]] - 
    Dq[p] Sin[eta] Sin[Theta[p]] Sinh[Zeta0[p]];
  y = Im[zcenter[p]] + (Cos[eta] Cosh[Zeta0[p]] Dq[p] Sin[Theta[p]] + 
      Cos[Theta[p]] Dq[p] Sin[eta] Sinh[Zeta0[p]]);
  
  f = (2.0/(Cosh[2*Zeta0[p]] - 
        Cos[2*eta])) (((Sin[eta]*Cosh[Zeta0[p]])^2)*
       means11[p, q] + ((Cos[eta]*Sinh[Zeta0[p]])^2)*means22[p, q] - 
      0.5*(Sin[2*eta]*Sinh[2*Zeta0[p]])*means12[p, q])]*)


(*In the global coordinate system*)
SEtaEtaAveG[p_, q_, zeta_, eta_] := Block[{f},
  x = Re[zcenter[p]] + Dq[p] Cos[eta] Cos[Theta[p]] Cosh[zeta] - 
    Dq[p] Sin[eta] Sin[Theta[p]] Sinh[zeta];
  y = Im[zcenter[p]] + (Cos[eta] Cosh[zeta] Dq[p] Sin[Theta[p]] + 
      Cos[Theta[p]] Dq[p] Sin[eta] Sinh[zeta]);
  
  f = (2.0/(Cosh[2*zeta] - Cos[2*eta])) (((Sin[eta]*Cosh[zeta])^2)*
       means11[p, q] + ((Cos[eta]*Sinh[zeta])^2)*means22[p, q] - 
      0.5*(Sin[2*eta]*Sinh[2*zeta])*means12[p, q])]


(*In the local yp coordinate system.*)
SEtaEtaAve[p_, q_, zeta_, eta_] := Block[{f},
  x = Dq[p] Cos[eta] Cosh[zeta];
  y = Dq[p] Sin[eta] Sinh[zeta];
  
  f = (2.0/(Cosh[2*zeta] - Cos[2*eta])) (((Sin[eta]*Cosh[zeta])^2)*
       means11[p, q] + ((Cos[eta]*Sinh[zeta])^2)*means22[p, q] - 
      0.5*(Sin[2*eta]*Sinh[2*zeta])*means12[p, q])]





(*Stored stress in the inclusion j due to the presence of crack with ttest number of micro-cracks.*)
StoredStressIncttt[j_,ttest1_]:=Block[{f,f1,f2,r,r1,rr0,rr1,rr2,rr},
Clear[Xi];
(*Print["I HAVE CHANGED {ii,ntot-ttest+1,ntot} TO {ii,ntot-ttest+1,ntot-1}, TO ELEMINATE THE EFFECT OF THE LAST MICROCRACK ON THE INCLUSION J."];*)
f  =  Table[If[ii!=j,ii,Unevaluated[Sequence[]]],{ii,Min[ntot-ttest1+1,ntot-1],Max[ntot-ttest1+1,ntot-1]}];
Print["coefficient test=  ", f];
f1 = Sum[means22[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];
f2 = Sum[means12[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];

rr0 = ((Dq[j] + r*Cos[eta0[j]])/l1[j])^2 + ((r*Sin[eta0[j]])/l2[j])^2;
  rr1 = r /. Solve[rr0 == 1, {r}];
  rr2 = If[rr1[[1]] > 0, rr1[[1]], rr1[[2]]];
  Print["r =  ", rr2];

r=Sqrt[2*Pi*(rr2)](Sigma22m[j]+I*Sigma12m[j]);
rr=SigmaEtaEta[j];
Xi=Zeta0[j]+I*eta0[j];
Print["SigmaEtaEta[j] =  ",rr];
Print["These Energy and SIF are not necessarily correct."];
Print["r = ",r];
r1 = f1 + I*f2 +r]


StoredStressInctt[j_,ttest1_]:=Block[{f,f1,f2,r,r1,rr0,rr1,rr2,rr},
Clear[Xi];
(*Print["I HAVE CHANGED {ii,ntot-ttest+1,ntot} TO {ii,ntot-ttest+1,ntot-1}, TO ELEMINATE THE EFFECT OF THE LAST MICROCRACK ON THE INCLUSION J."];*)
f  =  Table[If[ii!=j,ii,Unevaluated[Sequence[]]],{ii,Min[ntot-ttest1+1,ntot-1],Max[ntot-ttest1+1,ntot-1]}];
Print["coefficient test=  ", f];
f1 = Sum[means22[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];
f2 = Sum[means12[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];

(*rr0 = ((Dq[j] + r*Cos[eta0[j]])/l1[j])^2 + ((r*Sin[eta0[j]])/l2[j])^2;
  rr1 = r /. Solve[rr0 == 1, {r}];
  rr2 = If[rr1[[1]] > 0, rr1[[1]], rr1[[2]]];
  Print["r =  ", rr2];*)(*This rr0 is equal to dd[j](Cosh[Zeta0[j]+I*eta0temp]-1)*)

eta0temp=eta0[j];
r=Sqrt[2*Pi*(dd[j](Cosh[Zeta0[j]+I*eta0temp]-1))](Sigma22m[j]+I*Sigma12m[j]);
rr=SigmaEtaEta[j];
Xi=Zeta0[j]+I*eta0temp;
Print["SigmaEtaEta[j] =  ",rr];
Print["These Energy and SIF are not necessarily correct."];
Print["r = ",r];
r1 = f1 + I*f2 +r]


StoredStressInct[j_,ttest1_]:=Block[{f,f1,f2,r,r1,rr0,rr1,rr2,rr},
Clear[Xi];
(*Print["I HAVE CHANGED {ii,ntot-ttest+1,ntot} TO {ii,ntot-ttest+1,ntot-1}, TO ELEMINATE THE EFFECT OF THE LAST MICROCRACK ON THE INCLUSION J."];*)
f  =  Table[If[ii!=j,ii,Unevaluated[Sequence[]]],{ii,Min[ntot-ttest1+1,ntot-1],Max[ntot-ttest1+1,ntot-1]}];
Print["coefficient test=  ", f];
f1 = Sum[means22[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];
f2 = Sum[means12[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];

rr0 = ((Dq[j] + r*Cos[eta0[j]])/l1[j])^2 + ((r*Sin[eta0[j]])/l2[j])^2;
  rr1 = r /. Solve[rr0 == 1, {r}];
  rr2 = If[rr1[[1]] > 0, rr1[[1]], rr1[[2]]];
  Print["r =  ", rr2];

r=Sqrt[2*Pi*(rr2)]sigmmaLinemethodLength[j,e[j]];
rr=SigmaEtaEta[j];
Xi=Zeta0[j]+I*eta0[j];
Print["SigmaEtaEta[j] =  ",rr];
Print["These Energy and SIF are not necessarily correct."];
Print["r = ",r];
r1 = f1 + I*f2 +r]


eta0SIF[p_]:=Block[{f,r1,r2,k11,k22,q},If[ntot-ninitial!=ttest,q=p,q=p-ttest+1];
r1=k[p];
k11=Re[r1];
k22=Im[r1];
(*r1=ArcCos[(3*k22^2-Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];*)
r2=ArcCos[(3*k22^2+Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
f=If[k22>0,-r2,r2]]
(*For this one I have to write a new NewCenter[p_] one considering our new angle*)




eta0SIFt[p_]:=Block[{f,r1,r2,k11,k22,q},If[ntot-ninitial!=ttest,q=p,q=p-ttest+1];
r1=kcrack[p,q];
k11=Re[r1];
k22=Im[r1];
(*r1=ArcCos[(3*k22^2-Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];*)
r2=ArcCos[(3*k22^2+Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
f=If[k22>0,-r2,r2]]
(*For this one I have to write a new NewCenter[p_] one considering our new angle*)


End[] 
Protect@@Names["SIF`*"]
EndPackage[] 

BeginPackage["SIF`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","eta0p`","Addone`"}] 
EndPackage[] 











































