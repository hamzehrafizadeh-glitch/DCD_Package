(* ::Package:: *)

 BeginPackage[ "Fieldsxy`"];
 Unprotect@@Names["Fieldsxy`*"];
 ClearAll@@Names["Fieldsxy`*"];


SxxFinal::usage ="Sxx Stress Field.";
SyyFinal::usage ="Syy Stress Field.";
TxyFinal::usage ="Txy Stress Field.";
UxFinal::usage ="Ux Stress Field.";
UyFinal::usage ="Uy Stress Field.";
(*x::usage ="";
y::usage ="";*)



 Begin[ "Private`"];
 


(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["BCSolve.m"];
Needs["BCSolve`"];




W  = << W.dat;

Unprotect[x,y];

XY = << XY.dat;


dxdyXY[j_]:=Block[{f},f=Rotate2Dcc[(x-XY[[j,1]]),(y-XY[[j,2]]),XY[[j,4]]]]

dxXY[i_]:=Re[dxdyXY[i]]/.{Im[x]->0,Re[x]->x,Im[y]->0, Re[y]->y}
dyXY[i_]:=Im[dxdyXY[i]]/.{Im[x]->0,Re[x]->x,Im[y]->0, Re[y]->y}


asx=0;
asy=0;
asxx=0;
asyy=0;
atxy=0;
Clear[k,m]
k=Chi0;m=Mu0; 

For[i=1,i<NTOTAL+1,i++,
Clear[Ux1,Uy1];
Ux1=UxBE[dxXY[i],dyXY[i],W[[i]],W[[i+NTOTAL]],XY[[i,3]],k,m];Uy1=UyBE[dxXY[i],dyXY[i],W[[i]],W[[i+NTOTAL]],XY[[i,3]],k,m];
asx=asx+(Ux1*Cos[-XY[[i,4]]])+(Uy1*Sin[-XY[[i,4]]]);
asy=asy+(Uy1*Cos[-XY[[i,4]]])-(Ux1*Sin[-XY[[i,4]]]);]

For[i=1,i<NTOTAL+1,i++,
Clear[Sx1,Sy1,Txy1];
Sx1=SxBE[dxXY[i],dyXY[i],W[[i]],W[[i+NTOTAL]],XY[[i,3]],k];
Sy1=SyBE[dxXY[i],dyXY[i],W[[i]],W[[i+NTOTAL]],XY[[i,3]],k];
Txy1=TxyBE[dxXY[i],dyXY[i],W[[i]],W[[i+NTOTAL]],XY[[i,3]],k];

asyy=asyy+(Sx1*Sin[-XY[[i,4]]]^2)+(Sy1*Cos[-XY[[i,4]]]^2)-2*(Txy1*Cos[-XY[[i,4]]]*Sin[-XY[[i,4]]]);
asxx=asxx+(Sx1*Cos[-XY[[i,4]]]^2)+(Sy1*Sin[-XY[[i,4]]]^2)+2*(Txy1*Cos[-XY[[i,4]]]*Sin[-XY[[i,4]]]);
atxy=atxy-(Sx1-Sy1)*Cos[-XY[[i,4]]]*Sin[-XY[[i,4]]]+(Txy1*(Cos[-XY[[i,4]]]^2-Sin[-XY[[i,4]]]^2));]


SxxFinal=asxx;
SyyFinal=asyy;
TxyFinal=atxy;
UxFinal=asx;
UyFinal=asy;

SyyFinal >> SyyFinal.dat;
SxxFinal >> SxxFinal.dat;
TxyFinal >> TxyFinal.dat;

UxFinal >> UxFinal.dat;
UyFinal >> UyFinal.dat;



ClearSystemCache[];
 AbortKernels[];
CloseKernels[];


 End[];
 Protect@@Names["Fieldsxy`*"]
 EndPackage[]

BeginPackage["Fieldsxy`",{"GenerateMesh`","FieldsBE`","CoefMat`","BCSolve`"}] 
EndPackage[] 














































