(* ::Package:: *)

BeginPackage["PhiPsi`"]
Unprotect@@Names["PhiPsi`*"];
ClearAll@@Names["PhiPsi`*"];



Phi0::usage=" Phi is a scalar Muskhelishvili potential of complex variable method in plane Elasticity. zero suffix functions denote potentials in the matrix host."; 
Psi0::usage=" Phi and Psi are scalar Muskhelishvili potentials of complex variable method in plane Elasticity. zero suffix functions denote potentials in the matrix host."; 
Tau0::usage=" Traction field. zero suffix functions denote potentials in the matrix host."; 
TauFar::usage=" Traction field with far field source"; 
Alpha0::usage=" \[Alpha] = \[Sigma]11 - I \[Sigma]12, zero suffix functions denote potentials in the matrix host."; 
AlphaFar::usage=" \[Alpha] = \[Sigma]11 - I \[Sigma]12, Far suffix denote fields in the matrix host.";
PhiFarQ::usage=" Phi scalar potential associated with the far-field load at qth coordinate system."; 
PsiFarQ::usage=" Psi scalar potential associated with the far-field load at qth coordinate system.";
Beta0::usage=" \[Beta] = \[Sigma]11 + \[Sigma]22, zero suffix functions denote field with far field source"; 
BetaFar::usage=" \[Beta] = \[Sigma]11 + \[Sigma]22, zero suffix functions denote field with far field source"; 
(*v::usage=" v is a harmonic term, v = Exp[\[Xi]]"; *)
(*Xi::usage="  = \[Xi] = \[Zeta] + I\[Eta]"; 
XiQ::usage=" \[Xi]q[q] is a local q-th particle in the global z coordinate system. "; 
z::usage=" z= x + iy; coordinate of global system."; *)

Phi0z::usage=" Phi is a scalar Muskhelishvili potential of complex variable method in plane elasticity. zero suffix functions denote potentials in the matrix host. z suffix indicates that fields are calculated in z coordinate system."; 
Psi0z::usage=" Phi and Psi are scalar Muskhelishvili potentials of complex variable method in plane Elasticity. zero suffix functions denote potentials in the matrix host.z suffix indicates that fields are calculated in z coordinate system."; 
Alpha0z::usage=" \[Alpha] = \[Sigma]11 - I \[Sigma]12, zero suffix functions denote potentials in the matrix host. z suffix indicates that fields are calculated in z coordinate system."; 
AlphaFarz::usage=" \[Alpha] = \[Sigma]11 - I \[Sigma]12, Far suffix denote fields in the matrix host. z suffix indicates that fields are calculated in z coordinate system.";
PhiFarQz::usage=" Phi scalar potential associated with the far-field load at z coordinate system."; 
PsiFarQz::usage=" Psi scalar potential associated with the far-field load at z coordinate system.";
Beta0z::usage=" \[Beta] = \[Sigma]11 + \[Sigma]22, zero suffix functions denote field with far field source. z suffix indicates that fields are calculated in z coordinate system."; 
BetaFarz::usage=" \[Beta] = \[Sigma]11 + \[Sigma]22, zero suffix functions denote field with far field source. z suffix indicates that fields are calculated in z coordinate system."; 

PhiQInside::usage=" Phi is a scalar Muskhelishvili potential of complex variable method in plane elasticity. QInside suffix functions denote potentials in the inclusion. z suffix indicates that fields are calculated in z coordinate system."; 
PsiQInside::usage=" Phi and Psi are scalar Muskhelishvili potentials of complex variable method in plane Elasticity. QInside suffix functions denote potentials in the inclusion. z suffix indicates that fields are calculated in z coordinate system."; 
AlphaQ::usage=" \[Alpha] = \[Sigma]11 - I \[Sigma]12, zero suffix functions denote potentials in the matrix host. Q suffix indicates that fields are calculated in z coordinate system inside inclusion Q."; 
BetaQ::usage=" \[Beta] = \[Sigma]11 + \[Sigma]22, zero suffix functions denote field with far field source. Q suffix indicates that fields are calculated in z coordinate system inside inclusion Q."; 

Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
(* **Get["linearSolve.m"];*)
Needs["linearSolve`"];


(*PhiQs[q_] : 
  " Phi and Psi are scalar Muskhelishvili potentials of complex variable method in plane 
elasticity. Q suffix functions denote potentials calculated in "qth" coordinate system. 
s and r suffix denotes the singular and regular field calculations. PQr means calculating 
regular fields at qth coordinate system due to the presence of pth inclusion. FarQ suffix 
indicates potentials are generated due the far-field load at qth coordinate system."; 
PsiQs[q_]
PhiPQr[p_, q_]
PsiPQr[p_, q_]
PhiFarQ[q_]
PsiFarQ[q_]
Omegaprim[q_] : "It is derivative of Omega with respect to z."*)
(*(*v = Exp[Xi]*)
XiQ[q_] := 
 Module[{f, zq}, zq = z - zcenter[q]; f = ArcCosh[zq/dd[q]]]*)

PhiQs[q_] := Module[{f,v}, 
v = Exp[Xi]; f =  Sum[abmnpq`a[i, q] v^(-i) , {i, n}]] 

PsiQs[q_] := 
 Module[{f, f1, Psi0, Psi1,v}, 
   v = Exp[Xi];
  f1 = (Sinh[Zeta0[q]] /Sinh[Xi]) (v/v0[q] - 
      v0[q]/v); Psi0 = Sum[abmnpq`b[i, q] v^(-i), {i, n}]; 
  Psi1 = Sum[ i f1 abmnpq`a[i, q] v^(-i), {i, n}];
  f = Psi0 - Psi1] 

PhiPQr[p_, q_] := 
 Module[{f,v}, v = Exp[Xi]; f = Sum[anpq[i, p, q] (v^i + v^(-i)) , {i, n}]] 

PsiPQr[p_, q_] := Module[{f, Psi0 , Psi1, f1,v}, v = Exp[Xi];
  f1 = (Sinh[Zeta0[q]] Sinh[Xi - Zeta0[q]])/(Sinh[Xi]);
  Psi0 = 
   Sum[bnpqminus[i, p, q] v^(i) + bnpq[i, p, q] v^(-i) , {i, n}]; 
  Psi1 = Sum[2 i f1 anpq[i, p, q] ( v^(-i) - v^(i)), {i, n}] ;
  f =  Psi0 - Psi1 ]

PhiFarQ[q_] := 
 Module[{f,v}, v = Exp[Xi]; f =  Sum[faraa[q][[i]] (v^i + v^(-i)) , {i, n}]]

PsiFarQ[q_] := 
 Module[{f, f1, v}, v = Exp[Xi];
  f1 = (Sinh[Zeta0[q]]*Sinh[Xi - Zeta0[q]])/(Sinh[Xi]); 
  f = Sum[(farbb[q][[i]] - 2*i*f1*faraa[q][[i]])*v^(-i) 
     + (farbbminus[q][[i]] + 2*i*f1*faraa[q][[i]])*v^(i), {i, n}]] 

Omegaprim[q_] := 
 Module[{f}, f = Exp[-I Theta[q]] dd[q] Sinh[Xi]]

 Phi0[q_]:=Module[{f},f=PhiQs[q]+Sum[PhiPQr[p,q],{p,ntot}]]

Psi0[q_]:=Module[{f},f=PsiQs[q]+Sum[PsiPQr[p,q],{p,ntot}]]

(*I think I have to make conjugate of it for Tau formulas*)
Tau0[q_] := 
  Module[{f, f1, f2, Phizprim, Phizprimprim, Psizprim, zq, 
      yq  },
    f1 = Omegaprim[q] / Conjugate[Omegaprim[q]];
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phi0[q])\)\) / Omegaprim[q];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phizprim)\)\) / Omegaprim[q];
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Psi0[q])\)\)/Omegaprim[q];
     zq = dd[q] Cosh[Xi];
    yq = Exp[-I Theta[q]] zq; 
    f = 2.0 Mu0 (Phizprim + Conjugate[Phizprim] - 
            f1 (( Conjugate[yq] - 
                        
            yq) Phizprimprim - Phizprim + Psizprim) )]
TauFar[q_] := 
  Module[{f, f1, f2, Phizprim, Phizprimprim, Psizprim , 
   zq , 
      yq}, Clear[Psizprim];
    f1 = Omegaprim[q] / Conjugate[Omegaprim[q]];
    Phizprim = Simplify[ \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PhiFarQ[q])\)\) / Omegaprim[q]];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phizprim)\)\) / Omegaprim[q];
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PsiFarQ[q])\)\)/Omegaprim[q];
     zq = dd[q] Cosh[Xi]; 
    yq = Exp[-I Theta[q]] zq; 
    f = 2.0 Mu0 (Phizprim + Conjugate[Phizprim] - 
            f1 ((Conjugate[yq] - 
                        
            yq) Phizprimprim - Phizprim + Psizprim) )
    ]
Alpha0[q_] := 
  Module[{f, f1, f2, f3, Phizprim, Phizprimprim, Psizprim , 
      zq, yq },
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phi0[q])\)\) / Omegaprim[q];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phizprim)\)\) / Omegaprim[q];
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Psi0[q])\)\)/Omegaprim[q];
    f3 = Simplify[
        2.0 Mu0 (Phizprim + 
              Conjugate[Phizprim] - ((Conjugate[yq] - 
                        
            yq) Phizprimprim - Phizprim + Psizprim) )];
    zq = dd[q] Cosh[Xi]; 
    yq = Exp[-I Theta[q]] zq;
    f = f3
    ]

(*AlphaFar and TauFar are rong*)
AlphaFar0[q_] := 
  Module[{f, f1, f2, Phizprim, Phizprimprim, Psizprim, zq , yq },
    Clear[Xi];
    Phizprim = Simplify[ \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PhiFarQ[q])\)\)/ Omegaprim[q]];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phizprim)\)\) / Omegaprim[q];
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PsiFarQ[q])\)\)/Omegaprim[q];
     zq = dd[q]*Cosh[Xi]; 
    yq = Exp[-I*Theta[q]]*zq;
    f1 = 2.0*Mu0*(Phizprim + 
            Conjugate[Phizprim] - ((yq - Conjugate[yq]) Phizprimprim - Phizprim + Psizprim) ) 
    ]
AlphaFar[q_] := Module[{f}, f = FullSimplify[AlphaFar0[q]]]

Beta00[q_] := Module[{f1, f, Phizprim, zq}, 
    Clear[Xi];
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((Phi0[q])\)\) / Omegaprim[q];
    f = 4 Mu0 Re[(Phizprim + Conjugate[Phizprim])]]

Beta0[q_] := Module[{f}, f = Simplify[Beta00[q]]]

BetaFar0[q_] := Module[{f, Phizprim}, 
    Clear[Xi];
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(Xi\)]\((PhiFarQ[q])\)\) / Omegaprim[q];
    f = 4 Mu0 Re[(Phizprim + Conjugate[Phizprim])]]
BetaFar[q_] := Module[{f}, f = FullSimplify[BetaFar0[q]]]



















(*Calculating fields in the MATRIX MATERIAL and z coordinate system.*)
(*To calculate potentials and fields in the z coordinate system, we have to write 
potentials in the z coordinate system and use the singular terms. Singular potentials mean 
that we only consider the negative exponential terms. Here suffix z means that all 
calculations are converted into "z" coordinate system. Therefore, we can easily calculate 
the total fields in the z coordinate system.*)

PhiQsz[q_] := Module[{f,v}, 
v = Exp[XiQ[q]]; f =  Sum[abmnpq`a[i, q] v^(-i) , {i, n}]] 

PsiQsz[q_] := 
 Module[{f, f1, Psi0, Psi1,v}, 
   v = Exp[XiQ[q]];
  f1 = (Sinh[Zeta0[q]] /Sinh[XiQ[q]]) (v/v0[q] - 
      v0[q]/v); Psi0 = Sum[abmnpq`b[i, q] v^(-i), {i, n}]; 
  Psi1 = Sum[ i f1 abmnpq`a[i, q] v^(-i), {i, n}];
  f = Psi0 - Psi1] 

PhiFarQz[q_] := 
 Module[{f,v}, v = Exp[XiQ[q]]; f =  Sum[faraa[q][[i]] (v^i + v^(-i)) , {i, n}]]

PsiFarQz[q_] := 
 Module[{f, f1, v}, v = Exp[XiQ[q]];
  f1 = (Sinh[Zeta0[q]]*Sinh[XiQ[q] - Zeta0[q]])/(Sinh[XiQ[q]]); 
  f = Sum[(farbb[q][[i]] - 2*i*f1*faraa[q][[i]])*v^(-i) 
     + (farbbminus[q][[i]] + 2*i*f1*faraa[q][[i]])*v^(i), {i, n}]] 

(*Please note that here we just used the singular fields.*)
Phi0z:=Module[{f},f=Sum[PhiQsz[p],{p,ntot}]]

Psi0z:=Module[{f},f=Sum[PsiQsz[p],{p,ntot}]]

Alpha0z := 
  Module[{f, f1, f2, f3, Phizprim, Phizprimprim, Psizprim , 
      zq, yq },
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phi0z)\)\);
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phizprim)\)\) ;
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Psi0z)\)\);
    f3 = Simplify[
        2.0 Mu0 (Phizprim + 
              Conjugate[Phizprim] - ((Conjugate[z] - 
                        
            z) Phizprimprim - Phizprim + Psizprim) )];

    f = f3
    ]

AlphaFar0z[q_] := 
  Module[{f, f1, f2, Phizprim, Phizprimprim, Psizprim, zq , yq },
   
    Phizprim = Simplify[ \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PhiFarQz[q])\)\)];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phizprim)\)\) ;
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PsiFarQz[q])\)\);
    f1 = 2.0*Mu0*(Phizprim + 
            Conjugate[Phizprim] - ((z - Conjugate[z]) Phizprimprim - Phizprim + Psizprim) ) 
    ]
AlphaFarz[q_] := Module[{f}, f = FullSimplify[AlphaFar0z[q]]]

Alpha0z := 
  Module[{f, f1, f2, f3, Phizprim, Phizprimprim, Psizprim , 
      zq, yq },
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phi0z)\)\);
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phizprim)\)\) ;
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Psi0z)\)\);
    f3 = Simplify[
        2.0 Mu0 (Phizprim + 
              Conjugate[Phizprim] - ((Conjugate[z] - 
                        
            z) Phizprimprim - Phizprim + Psizprim) )];

    f = f3
    ]

AlphaFar0z[q_] := 
  Module[{f, f1, f2, Phizprim, Phizprimprim, Psizprim, zq , yq },
   
    Phizprim = Simplify[ \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PhiFarQz[q])\)\)];
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phizprim)\)\) ;
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PsiFarQz[q])\)\);
    f1 = 2.0*Mu0*(Phizprim + 
            Conjugate[Phizprim] - ((z - Conjugate[z]) Phizprimprim - Phizprim + Psizprim) ) 
    ]
AlphaFarz[q_] := Module[{f}, f = FullSimplify[AlphaFar0z[q]]]

Beta00z := Module[{f1, f, Phizprim, zq}, 
    
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phi0z)\)\) ;
    f = 4 Mu0 Re[(Phizprim + Conjugate[Phizprim])]]

Beta0z := Module[{f}, f = Simplify[Beta00z]]

BetaFar0z[q_] := Module[{f, Phizprim}, 
   
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PhiFarQz[q])\)\);
    f = 4 Mu0 Re[(Phizprim + Conjugate[Phizprim])]]
BetaFarz[q_] := Module[{f}, f = FullSimplify[BetaFar0z[q]]]




(*Calculating fields in the INCLUSIONS and z coordinate system.*)
(*To calculate potentials and fields in the z coordinate system, we have to write 
potentials in the z coordinate system and using c and d expansion coefficients. Therefore, we can easily calculate 
the total fields in the z coordinate system.*)

Clear[PhiQInside,PsiQInside,AlphaQ,BetaQ0]
(*dminus[i_,q_]:=Module[{f},f=abmnpq`d[i,q]-2 i Sinh[2*Zeta0[q]] abmnpq`c[i,q]]*)
PhiQInside[q_] := Module[{f,v}, 
v = Exp[XiQ[q]]; f =  Sum[abmnpq`c[i, q] (v^i + v^(-i)) , {i, n}]] 

PsiQInside[q_] := 
 Module[{f, f1, Psi0, Psi1,v}, 
   v = Exp[XiQ[q]];
  f1 = (Sinh[Zeta0[q]] /Sinh[XiQ[q]]) (v/v0[q] - v0[q]/v);Psi0 = Sum[abmnpq`dminus[i, q] v^(i) + abmnpq`d[i, q] v^(-i), {i, n}]; 
  Psi1 = Sum[ i f1 abmnpq`c[i, q] ( v^(-i) - v^(i)), {i, n}];
  f = Psi0 - Psi1] 

AlphaQ[q_] := 
  Module[{f, f1, f2, f3, Phizprim, Phizprimprim, Psizprim , 
      zq, yq },
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PhiQInside[q])\)\) ;
    Phizprimprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((Phizprim)\)\) ;
    Psizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PsiQInside[q])\)\);
    f3 = Simplify[
        2.0 MuQ[q] (Phizprim + 
              Conjugate[Phizprim] - ((Conjugate[z] - 
                        
            z) Phizprimprim - Phizprim + Psizprim) )];
 
    f = f3
    ]
BetaQ0[q_] := Module[{f1, f, Phizprim, zq}, 
    Phizprim = \!\(
\*SubscriptBox[\(\[PartialD]\), \(z\)]\((PhiQInside[q])\)\) ;
    f = 4 MuQ[q] Re[(Phizprim + Conjugate[Phizprim])]]

BetaQ[q_] := Module[{f}, f = Simplify[BetaQ0[q]]]







(* ::Input:: *)
(**)


End[] 
Protect@@Names["PhiPsi`*"]
EndPackage[] 

BeginPackage["PhiPsi`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" ,"TemperatureInclusions`", "farfield`" ,"linearSolve`"}] 
EndPackage[] 




