(* ::Package:: *)

BeginPackage["farfield`"]
Unprotect @@ Names["farfield`*"];
ClearAll @@ Names["farfield`*"];

aa::usage=" aa"; 
bbm::usage=" bbm"; 
faraa::usage=" faraa"; 
farbb::usage=" faraa"; 
farbbminus::usage=" farbbminus"; 
faraaminus::usage=" faraaminus"; 
bcc::usage=" bcc"; 
bcbar::usage=" bcbar"; 

Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["matrix.m"];
Needs["matrix`"];



ClearAll[listaa,listaaminus,listbb,listbbminus]


aa[q_]:=Module[{f}, f=Exp[-I Theta[q]](dd[q]/(16 Mu0))(s11+s22)]

bbm[q_] :=Module[{f}, f= aa[q]v0[1]^(-2) +  Exp[I Theta[q]](dd[q]/(8 Mu0))(s22-s11 + 2 I s12)]

(*"faraa"*)
faraa[q_]  := Module[{f}, f= Table[If[i ==1&&q<= ntot, Exp[-I Theta[q]](dd[q]/(16  Mu0))(s11+s22),0],{i,1,n}]]

(*"farbbminus"*)
farbbminus[q_]:= Module[{f, f1, f2, f3,f4}, 
f1=Exp[-I Theta[q]](dd[q]/(16*Mu0))(s11+s22);
f2= f1 v0[q]^(-2) +  Exp[I Theta[q]](dd[q]/(8*Mu0))(s22-s11 + 2 I s12) ;
f= Table[If[i ==1&&q<= ntot, f2,0],{i,1,n}]]



(*"faraaminus"*)
r5[q_]  := Module[{f}, f= faraa[q]]
faraaminus[q_]:=Module[{f},f=r5[q]]

(*"farbb"*)
farbb[q_] := Module[{f}, f= Table[ farbbminus[q][[i]] + 2 i Sinh[2*Zeta0[q]] faraa[q][[i]],{i,1,n}]]

listaa[q_]  := Module[{f},f= Table[faraa[q][[i]] , {i,n}]]
 listaaminus[q_] := Module[{f},f=listaa[q]]

listbbminus[q_]  := Module[{f},f= Table[farbbminus[q][[i]] , {i,n}]]
listbb[q_] := Module[{f},f= Table[farbb[q][[i]], {i,n}]]
(*"-----------VECTOR OF BOUNDARY CONDITIONS--------------"*)


ClearAll[bcf1, bcf2, bcf3, bcf,bc]

bcf1[i_, j_,q_] := Module[{f},If[i == 1,f= -Chi0 listaa[q][[j]] + v0[q]^(2 j)Conjugate[listbbminus[q][[j]]],0]
]

bcf2[i_, j_,q_] := Module[{f},If[i == 2,f =Chi0 v0[q]^(2j) listaaminus[q][[j]]- Conjugate[listbb[q][[j]]],0]
]

bcf3[i_, j_,q_] := Module[{f},If[i == 3,f =-listaa[q][[j]] - v0[q]^(2 j) Conjugate[listbbminus[q][[j]]],0]
]

bcf4[i_, j_,q_] := Module[{f},If[i == 4,f =- v0[q]^(2 j) listaaminus[q][[j]] - Conjugate[listbb[q][[j]]],0]
]

bcf[i_, j_,q_]:= Module[{f}, f =( bcf1[i,j,q]+ bcf2[i,j,q]+bcf3[i,j,q]+bcf4[i,j,q])]

ClearAll[bbc]

bbc[j_,q_]:=Module[{f}, f=Table[bcf[j,i,q], {k,1},{i,n}]]
Clear[bc]


bc[q_] :=Module[{f}, f=Table[bbc[i,q],{i,4}]]


bbctest[j_,q_]:=Module[{f}, f=Table[bcf[i,j,q], {k,1},{i,4}]]
bcc :=Module[{f}, f=Flatten[Table[bbctest[i,q],{q,ntot},{i,n}]]]

bcr=Simplify[(Re[bcc])];
bci=Simplify[Im[bcc]];
bcbar=Flatten[{{bcr},{bci}}];

End[] 
Protect@@Names["farfield`*"]
EndPackage[] 

BeginPackage["farfield`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`"}] 
EndPackage[] 









