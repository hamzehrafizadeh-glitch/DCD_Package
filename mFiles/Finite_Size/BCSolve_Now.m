(* ::Package:: *)

 BeginPackage[ "BCSolve`"];
 Unprotect@@Names["BCSolve`*"];
 ClearAll@@Names["BCSolve`*"];


Sij::usage ="Sxx Stress Field.";
Sbc::usage ="Syy Stress Field.";



 Begin[ "Private`"];
 


AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]]
Get["CoefMat.m"];
Needs["CoefMat`"];
XY = << XY.dat;
Cij = << Cij.dat;
Uij = << Uij.dat;
SMM = <<SMM.dat;



(*If we have any mixed or pure displacement boundary condition,
the Sij matrix should be changed according to the boundary condition.*)

(*For axample*)
(*n1a=2;
n1=5;
Sijt=Table[0,{i,1,2*NTOTALt},{j,1,2*NTOTALt}];
"3a"
Do[Sijt[[i,j]]=Cijt[[i,j]],{i,1,2*NTOTALt},{j,1,2*NTOTALt}];
"3b"
Do[Sijt[[i ,j]]=Uijt[[i+NTOTALt,j]],{i,n1a+1,n1},{j,1,2*NTOTALt}];
Do[Sijt[[i+NTOTALt,j]]=Uijt[[i,j]],{i,n1a+1,n1},{j,1,2*NTOTALt}];
(******if both XY are fixed******)"3c"*)


(*(*Pure tension*)
(*------>>>   Form Solutions MATRIX    <<<---~-----------*)


Sij = Table[0, {i, 1, 2*NTOTAL}, {j, 1, 2*NTOTAL}];
Do[Sij[[i, j]] = Cij[[i, j]], {i, 1, 2*NTOTAL}, {j, 1, 2*NTOTAL}]

(*----------------------External Loading-----------------*)

S = Table[0, {i, 2*NTOTAL}];
SFar = Table[0, {i, 2*NTOTAL}];
(*Set1*)
Do[SFar[[i]] = 1, {i, 1, n1}];
Do[SFar[[i]] = 0, {i, NTOTAL + 1, NTOTAL+n1}];

(*set2*)
Do[SFar[[i]] = 1, {i, n1 + 1,2n1}];
Do[SFar[[i]] = 0, {i, NTOTAL +n1 + 1, NTOTAL +2n1}];

(*set3*)
Do[SFar[[i]] = 0, {i, 2n1+n2+1,NTOTAL-4*ncorner}];
Do[SFar[[i]] = 0, {i,NTOTAL + 2n1+n2+1,2*NTOTAL-4*ncorner}];

(*set4*)

Do[SFar[[i]] = 0, {i, 2n1+1,2n1+n2}];
Do[SFar[[i]] = 0, {i,NTOTAL + 2n1+2,NTOTAL+2n1+n2}];*)




(*------>>>   Form Solutions MATRIX    <<<---~-----------*)

(*Tension and 2 fixed displacement points*)
Sij = Table[0, {i, 1, 2*NTOTAL}, {j, 1, 2*NTOTAL}];
Do[Sij[[i, j]] = Cij[[i, j]], {i, 1, 2*NTOTAL}, {j, 1, 2*NTOTAL}]

(*Two fixed points are:
Table[{XY[[i,1]],XY[[i,2]]},{i,{2n1+n2+1,2n1+2n2}}]*)
Do[Sij[[i (*+NTOTAL ****** if both XY are fixed ******* *),j]]=Uij[[i,j]],{i,{2n1+n2+1,2n1+2n2}},{j,1,2*NTOTAL}];
Do[Sij[[i+NTOTAL,j]]=Uij[[i+NTOTAL,j]],{i,{(*2n1+n2+1,*)2n1+2n2}},{j,1,2*NTOTAL}]; (*   ***** if both XY are fixed ***** *)  

(*Do[Sij[[i (*+NTOTAL ****** if both XY are fixed ******* *),j]]=Uij[[i,j]],{i,{1,n1+1}},{j,1,2*NTOTAL}];
Do[Sij[[i+NTOTAL,j]]=Uij[i+NTOTAL,j],{i,{1,n1+1}},{j,1,2*NTOTAL}]; (*   ***** if both XY are fixed ***** *)  *)

(*Do[Sij[[i+NTOTAL,j]]=Uij[i+NTOTAL,j],{i,{1}},{j,1,2*NTOTAL}]; *)


(*----------------------External Loading-----------------*)

S = Table[0, {i, 2*NTOTAL}];
SFar = Table[0, {i, 2*NTOTAL}];

(*I meight better to write 
SFar1=Table[t=-XY[[j,4]]+Pi/2;Tnfield[t],{j,NTOTAL}];

SFar2=Table[t=-XY[[j,4]]+Pi/2;Tsfield[t],{j,NTOTAL}];
SFar=Join[SFar1,SFar2];*)


(*Set1*)
(*Normal*)Do[SFar[[i]] = 1, {i, 1, n1}];
(*Traction*)Do[SFar[[i+NTOTAL]] = 0, {i, 1,n1}];

(*set2*)
(*Normal*)Do[SFar[[i]] = 1, {i, n1 + 1,2n1}];
(*Traction*)Do[SFar[[i+NTOTAL]] = 0, {i, n1 + 1, 2n1}];

(*set3*)
(*Normal*)Do[SFar[[i]] = 0, {i, 2n1+n2+2,NTOTAL-4*ncorner-1}];
(*Traction*)Do[SFar[[i+NTOTAL]] = 0, {i, 2n1+n2+1,NTOTAL-4*ncorner}];

(*set4*)

(*Normal*)Do[SFar[[i]] = 0, {i, 2n1+1,2n1+n2}];
(*Traction*)Do[SFar[[i+NTOTAL]] = 0, {i, 2n1+2,2n1+n2}];

(*Fixed Displacements.*)
Do[SFar[[i ]]=0,{i,{2n1+n2+1,2n1+2n2}}];
Do[SFar[[i+NTOTAL]]=0,{i,{2n1+n2+1,2n1+2n2}}]; (*   ***** if both XY are fixed ***** *) 

Do[SFar[[i ]]=1,{i,{(*1,*)n1+1}}]; (*For the case of soft Incusion, I have put this one zero.*)
SFar[[1 ]]=0; (*For the case of hard Incusion, I have put this one zero.*)
Do[SFar[[i+NTOTAL]]=0,{i,{1,n1+1}}]; (*   ***** if both XY are fixed ***** *)  

SFar[[NTOTAL]]=0.5
SFar[[NTOTAL-1]]=0.5  (*Displacement*)
SFar[[NTOTAL-2]]=0.5
SFar[[NTOTAL-3]]=0.5

SFar[[2*NTOTAL]]= 0.5
SFar[[2NTOTAL-1]]=-0.5  (*Displacement*)
SFar[[2NTOTAL-2]]=0.5;0  (*For the case of soft Incusion, I have put this one zero.*)
SFar[[2NTOTAL-3]]=0;-0.5 (*For the case of hard Incusion, I have put this one zero.*)


Print["SFar",SFar];





(*S=SFar;*)
S=SFar-SMM;
S>>S.dat;
(*----Solving the system------*)

 Solutions = LinearSolve[Sij, S];
Print["Solved"];
 Solutions >> W.dat;


(*Set1
Clear[s22,sx1,sx2,txy1,sqw,s12]
sx1=0;sy1=1;sqw=0;txy1=0;
s22=(sx1*Sin[sqw]^2+sy1*Cos[sqw]^2-2*txy1*Sin[sqw]*Cos[sqw])
s12=-(sx1-sy1)*Sin[sqw]*Cos[sqw]+txy1*(Cos[sqw]^2-Sin[sqw]^2)

set2
Clear[s22,sx1,sx2,txy1,sqw,s12]
sx1=0;sy1=1;sqw=-Pi;txy1=0;
s22=(sx1*Sin[sqw]^2+sy1*Cos[sqw]^2-2*txy1*Sin[sqw]*Cos[sqw])
s12=-(sx1-sy1)*Sin[sqw]*Cos[sqw]+txy1*(Cos[sqw]^2-Sin[sqw]^2)

set3
Clear[s22,sx1,sx2,txy1,sqw,s12]
sx1=0;sy1=0;sqw=-Pi/2;txy1=1;
s22=(sx1*Sin[sqw]^2+sy1*Cos[sqw]^2-2*txy1*Sin[sqw]*Cos[sqw])
s12=-(sx1-sy1)*Sin[sqw]*Cos[sqw]+txy1*(Cos[sqw]^2-Sin[sqw]^2)


set4
Clear[s22,sx1,sx2,txy1,sqw,s12]
sx1=0;sy1=0;sqw=Pi/2;txy1=1;
s22=(sx1*Sin[sqw]^2+sy1*Cos[sqw]^2-2*txy1*Sin[sqw]*Cos[sqw])
s12=-(sx1-sy1)*Sin[sqw]*Cos[sqw]+txy1*(Cos[sqw]^2-Sin[sqw]^2)*)


 End[];
 Protect@@Names["BCSolve`*"]
 EndPackage[]

BeginPackage["BCSolve`",{"GenerateMesh`","FieldsBE`","CoefMat`"}] 
EndPackage[] 




















































