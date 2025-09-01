(* ::Package:: *)

BeginPackage["abmnpq`"]
Unprotect @@ Names["abmnpq`*"];
ClearAll @@ Names["abmnpq`*"];


Eta::usage=" Eta"; 
Mu::usage="Mu"; 

amp::usage=" A_mp = Series coefficients; m determines the harmonic terms (\!\(\*SuperscriptBox[\(v\), \(-n\)]\)) and p, the inclusion. "; 
bmp::usage=" B_mp = Series coefficients; m determines the harmonic terms (\!\(\*SuperscriptBox[\(v\), \(n\)]\)) and p, the inclusion. "; 
cmp::usage=" C_mp = Series coefficients; m determines the harmonic terms (\!\(\*SuperscriptBox[\(v\), \(-n\)]\)) and p, the inclusion. "; 
dmp::usage=" D_mp = Series coefficients; m determines the harmonic terms (\!\(\*SuperscriptBox[\(v\), \(n\)]\)) and p, the inclusion. "; 
a::usage=" List of A_mp components"; 
b::usage=" List of B_mp components"; 
c::usage=" List of C_mp components"; 
d::usage=" List of D_mp components"; 

smnpq::usage=" smnpq[m,n,p,q]"; 
bmnpqa::usage=" bmnpqa[m,n,p,q] To have a better understanding of the coefficient amnpq and bmnpq, please refer to ESA report. a and b reffers to the coefficients of the harmonic \!\(\*SuperscriptBox[\(v\), \(-n\)]\) and \!\(\*SuperscriptBox[\(v\), \(n\)]\) terms respectively. "; 
bmnpqb::usage=" bmnpqb[m,n,p,q]To have a better understanding of the coefficient amnpq and bmnpq, please refer to ESA report. a and b reffers to the coefficients of the harmonic \!\(\*SuperscriptBox[\(v\), \(-n\)]\) and \!\(\*SuperscriptBox[\(v\), \(n\)]\) terms respectively. "; 
anpq::usage=" anpq[m,n,p,q] To have a better understanding of the coefficient amnpq and bmnpq, please refer to ESA report. a and b reffers to the coefficients of the harmonic \!\(\*SuperscriptBox[\(v\), \(-n\)]\) and \!\(\*SuperscriptBox[\(v\), \(n\)]\) terms respectively. "; 
bnpq::usage=" bnpq[m,n,p,q] To have a better understanding of the coefficient amnpq and bmnpq, please refer to ESA report. a and b reffers to the coefficients of the harmonic \!\(\*SuperscriptBox[\(v\), \(-n\)]\) and \!\(\*SuperscriptBox[\(v\), \(n\)]\) terms respectively. "; 
bnpqminus::usage=" bnpqminus[m,n,p,q] To have a better understanding of the coefficient amnpq and bmnpq, please refer to ESA report. a and b reffers to the coefficients of the harmonic \!\(\*SuperscriptBox[\(v\), \(-n\)]\) and \!\(\*SuperscriptBox[\(v\), \(n\)]\) terms respectively. "; 

Begin["`Private`"] 

(*Change the path into mFiles path.*)
(*AppendTo[$Path,FileNameJoin[{FilePath, "mFiles/Multipole_Method"}]];*)

(* ***Get["expansionCoefficients.m"];*)
Needs["ExpansionCoefficient`"];


Eta = <<Eta.1;
Mu = <<Mu.1;

amp:= Module[{f},f=Table[a[i,j],{i,n},{j,ntot }]]
bmp:= Module[{f},f=Table[b[i,j],{i,n},{j,ntot }]]
cmp:= Module[{f},f=Table[c[i,j],{i,n},{j,ntot }]]
dmp:= Module[{f},f=Table[d[i,j],{i,n},{j,ntot }]]

smnpq[m_,i_,p_,q_]:=Block[{f},f= Eta[[m,Abs[i],p,q]]Exp[-I Theta[p,q]]]
bmnpqa[m_,i_,p_,q_]:=Block[{f}, f=If[p==q,0,(-(dd[p]/2) (m/(m+1)) ((v0[p]-(1/v0[p]))^2) Exp[I Theta[p,q]] Mu[[m+1,Abs[i],p,q]]+(Abs[i] Exp[-I Theta[p,q]] (-v0[q]^(-2)+1)-m Exp[I Theta[p,q]] (-v0[p]^(-2)+1)) Eta[[m,Abs[i],p,q]]+((v0[q]-(1/v0[q]))^2) Exp[-I Theta[p,q]] (Sum[(2 k+Abs[i]) Eta[[m,2 k+Abs[i],p,q]],{k,0,nmax}])+(Exp[I Theta[q]] Conjugate[zpq[p,q]]-zpq[p,q] Exp[-I Theta[q]]) Exp[I Theta[p]] Mu[[m,Abs[i],p,q]]+((Mu[[m+1,Abs[i],p,q]]+If[m==1,MuClosePQ[m-1,Abs[i],p,q],Mu[[m-1,Abs[i],p,q]]]) (dd[p]/2)-Eta[[m,Abs[i],p,q]]) (Exp[-I Theta[p,q]]-Exp[I Theta[p,q]]))]]
bmnpqb[m_,i_,p_,q_]:=Block[{f}, f=Eta[[m,Abs[i],p,q]] Exp[I*Theta[p,q]]]
anpq[i_,p_,q_]:=Module[{f,f1,f2},f1=Sum[amp[[m,p]] smnpq[m,i,p,q],{m,n}]]
bnpq[i_,p_,q_]:=Block[{f},f=If[p==q,0,Sum[bmp[[m,p]]*bmnpqb[m,i,p,q]+amp[[m,p]]*bmnpqa[m,i,p,q],{m,n}]]]
bnpqminus[i_,p_,q_]:=Module[{f},f=bnpq[i,p,q]-2 i Sinh[2*Zeta0[q]] anpq[i,p,q]]


End[] 
Protect @@ Names["abmnpq`*"]
EndPackage[] 

BeginPackage["abmnpq`",{"geometry`", "ExpansionCoefficient`"}] 
EndPackage[] 










































