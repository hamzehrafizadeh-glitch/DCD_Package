(* ::Package:: *)

BeginPackage["DisplacementFields`"]
Unprotect@@Names["DisplacementFields`*"];
ClearAll@@Names["DisplacementFields`*"];


ux::usage=" Displacement field in zq local coordinate system"; 
uy::usage=" Displacement field in zq local coordinate system"; 
uxz::usage=" Displacement field in z global coordinate system"; 
uyz::usage=" Displacementfield in z global coordinate system"; 
ufarz::usage="ufarz[q_]";
u0xz::usage="ufarz[q_]";
u0yz::usage="ufarz[q_]";
uQz::usage="uQz[q_]: Displacement inside the inclusion";

uxfarz::usage=" Displacement field in z global coordinate system"; 
uyfarz::usage=" Displacementfield in z global coordinate system"; 
(*x::usage=" "; 
y::usage=" "; *)


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)

(* ***Get["StressFields.m"];*)
Needs["StressFields`"];
Needs["CombinDCDTEMP`"];


Unprotect[z,Xi,x,y,ux,uy,uxz,uyz];
Clear[ux,uy,uxz,uyz]

Omegaprim[q_] := 
 Module[{f}, f = Exp[-I Theta[q]] dd[q] Sinh[Xi]]

(*Phi0 and Psi0 are scalar potentials in the matrix material.*)

u0z[q_] := Module[{f, phizprim,zq,yq},
 Clear[Xi];
phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phi0[q])\)\) / Omegaprim[q];
 zq = dd[q] Cosh[Xi]; 
      yq = Exp[-I Theta[q]] zq;
f=(Chi0*Phi0[q] - (yq - Conjugate[yq] + zcenter[q] - 
            Conjugate[
             zcenter[q]])(Conjugate[phizprim]) - Conjugate[Psi0[q]] )]


ufarz[q_] := Module[{f, phizprim,zq,yq},
phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PhiFarQ[q])\)\) / Omegaprim[q];
zq = dd[q] Cosh[Xi]; 
      yq = Exp[-I Theta[q]] zq;
f=(Chi0 PhiFarQ[q] - (yq - Conjugate[yq] + zcenter[q] - 
            Conjugate[
             zcenter[q]])(Conjugate[phizprim]) - Conjugate[PsiFarQ[q]] )]
uout[q_]:=ufarz[q]+u0z[q]

ux[q_] := Module[{f}, f = Re[  uout[q]]]
uy[q_] := Module[{f}, f = Im[uout[q]]]


(*phiQ and psiQ are scalar potentials inside the inclusion q.*)
psiQ[q_] := Module[{f, \[Psi]0 , \[Psi]1,v}, 
v = Exp[Xi];
\[Psi]0 = Sum[(abmnpq`d[i, q]- 2 i Sinh[2 Zeta0[q]] abmnpq`c[i, q])v^(i)+abmnpq`d[i, q]v^(-i) ,{i,n}]; 
\[Psi]1=-2 Sum[(i abmnpq`c[i, q]( v^(i)-v^(-i))),{i, n}]((Sinh[Zeta0[q]] Sinh[Xi-Zeta0[q]])/(Sinh[Xi])) ;
f =  \[Psi]0 - \[Psi]1 ]

phiQ[q_] := Module[{f,v}, 
v = Exp[Xi];f=  Sum[abmnpq`c[i, q]( v^(i) + v^(-i)),{i, n}] ] 
Print["\[Phi]1"]

(*uQz[q_]: Displacement inside the inclusion*)
uQz[q_] := Module[{f, phizprim,zq,yq},
 Clear[Xi];
phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((phiQ[q])\)\) / Omegaprimt[q];
 zq = dd[q] Cosh[Xi]; 
      yq = Exp[-I Theta[q]] zq;
f=(Chi[q]*phiQ[q] - (yq - Conjugate[yq] + zcenter[q] - 
            Conjugate[
             zcenter[q]])(Conjugate[phizprim]) - Conjugate[psiQ[q]] )]
             




uxz[i_] := 
 Module[{f, f1, f2, f3}, Clear[x, y, zp, zq, Xi, z];
  f1 =ux[i];
  Xi =XiQ[i];
  f2 = f1;
z = x + I y;
  f3 = f2
  ]

uyz[i_] := 
 Module[{f, f1, f2, f3}, Clear[x, y, zp, zq, Xi, z];
  f1 =uy[i];
  Xi = XiQ[i];
  f2 = f1;
z = x + I y;
  f3 = f2
  ]


u0xz[i_] := 
 Module[{f, f1, f2, f3}, Clear[x, y, zp, zq, Xi, z];
  f1 =Re[u0z[i]];
  Xi =XiQ[i];
  f2 = f1;
z = x + I y;
  f3 = f2
  ]

u0yz[i_] := 
 Module[{f, f1, f2, f3}, Clear[x, y, zp, zq, Xi, z];
  f1 =Im[u0z[i]];
  Xi = XiQ[i];
  f2 = f1;
z = x + I y;
  f3 = f2
  ]



End[] 
Protect@@Names["DisplacementFields`*"]
EndPackage[] 

(*BeginPackage["DisplacementFields`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","CombinMMBE`","CombinDCDTEMP`"}] 
EndPackage[] *)

BeginPackage["DisplacementFields`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","CombinDCDTEMP`"}] 
EndPackage[] 




















