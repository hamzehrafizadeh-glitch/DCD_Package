(* ::Package:: *)

BeginPackage["Temperature1Inclusion`"]
Unprotect@@Names["Temperature1Inclusion`*"];
ClearAll@@Names["Temperature1Inclusion`*"];



(*In this file thermal stress fields are calculated for circular shape inclusions. 

Then fields in z coordinate system are transformed to the local zq coordinate system 
and then it is rotated to the yq local coordinate system 
(SigmaxyYq, SigmaxxYq, SigmayyYq). 

SEtaEta: is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function.*)






Bettal0 ::usage ="For Matrix material: 2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain; 2.*ThermalExpansion0; is for plane stress";
BettalQ::usage ="For Inclusions: 2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain; 2.*ThermalExpansion0; is for plane stress";
Fn::usage ="It is a series of the Expansion coefficient of the Far-field complex temperature function.";
Resnq::usage ="Resnq[n_,q_]: The real part of the expansion coefficient of the disturbance temperature field induced by the inclusion in the matrix material vanishes at infinity. ";
Imsnq::usage ="Imsnq[n_,q_]: The imaginary part of the expansion coefficient of the disturbance temperature field induced by the inclusion in the matrix material vanishes at infinity. ";
ReGnq::usage ="ReGnq[n_,q_]: The real part of the expansion coefficient of the disturbance temperature field in the inclusion. ";
ImGnq::usage ="ImGnq[n_,q_]: The imaginary part of the expansion coefficient of the disturbance temperature field in the inclusion. ";
OmegaFar::usage ="OmegaFar[q_]: Far-field complex temperature function.";
OmegaSQ::usage ="OmegaSQ[q_]: The disturbance temperature field induced by the inclusion in the matrix material and it is vanishing at infinity.";
OmegaGQ::usage ="OmegaGQ[q_]: The complex temperature function inside the inclusion.";
TFar::usage ="Far-field temperature function";
TMatrix::usage ="Temperature field in the matrix material.";
TQ::usage ="Temperature field in the inclusion.";
FnFinal::usage ="Expansion coefficients of the Far-field temperature. ";
SnqFinal::usage ="Expansion coefficients of the disturbance temperature field induced by the inclusion in the matrix material vanishes at infinity. ";
GnqFinal::usage ="Expansion coefficients of the disturbance temperature field induced by the inclusion inside the inclusion. ";
TempratureInteraction1::usage ="It is an additional term in the boundary condition (1) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";
TempratureInteraction2::usage ="It is an additional term in the boundary condition (2) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";


Begin["`Private`"] ;





(*If you need to run the file seperatly, you need to set the path*)
(**AppendTo[$Path,FileNameJoin[{FilePath, "mFiles/Multipole_Method"}]];**)
(*Print["FilePath = ",$Path]*)
Get["geometry.m"];
Needs["geometry`"]
Print["Temperature_1Inclusion = ",TS];


Bettal0 = 2.*ThermalExpansion0*(1 + Nu0);2.*ThermalExpansion0;

BettalQ[q_]:=Module[{f,f1},f1=2.*ThermalExpansionQ[q];f=2.*ThermalExpansionQ[q]*(1 + Nu[q])]

(*2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain
 2.*ThermalExpansion0; is for plane stress *)


gamma=5.+I*2;
Fn:={dd[q]*Conjugate[gamma]/2.}
ReFn[n_]:= Module[{f}, f= Re[Fn[[n]]]];
ImFn[n_]:= Module[{f}, f= Im[Fn[[n]]]];



Resnq[n_,q_] :=Module[{f,f1,f2},
f1=((ThermalConductivityQBar[q]-1)*(1+v0[q]^(2n)));
f2= (ThermalConductivityQBar[q]+(v0[q]^(n)+v0[q]^(-n))/(v0[q]^(n)-v0[q]^(-n)));
f= -ReFn[n]*f1/f2];
Imsnq[n_,q_] :=Module[{f,f1,f2},
f1=((ThermalConductivityQBar[q]-1)*(1-v0[q]^(2n)));
f2= (ThermalConductivityQBar[q]+(v0[q]^(n)-v0[q]^(-n))/(v0[q]^(n)+v0[q]^(-n)));
f= -ImFn[n]*f1/f2]

ReGnq[n_,q_] :=Module[{f,f1,f2,f3},
f1=((ThermalConductivityQBar[q]-1)*(1+v0[q]^(2n)));
f2= (ThermalConductivityQBar[q]+(v0[q]^(n)+v0[q]^(-n))/(v0[q]^(n)-v0[q]^(-n)));
f3=(1+v0[q]^(-2n))/(v0[q]^(2n)-v0[q]^(-2n));
f= ReFn[n]*(1+f1*f3/f2)/ThermalConductivityQBar[q]];
ImGnq[n_,q_] :=Module[{f,f1,f2,f3},
f1=((ThermalConductivityQBar[q]-1)*(1-v0[q]^(2n)));
f2= (ThermalConductivityQBar[q]+(v0[q]^(n)-v0[q]^(-n))/(v0[q]^(n)+v0[q]^(-n)));
f3=(-1+v0[q]^(-2n))/(v0[q]^(2n)-v0[q]^(-2n));
f= ImFn[n]*(1+f1*f3/f2)/ThermalConductivityQBar[q]]



(*Complex temperature function*)
OmegaFar[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
f1=Length[Fn];
f=Sum[Fn[[i]](v^(i)+v^(-i)),{i,f1}]]

OmegaSQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
f1=Length[Fn];
f=Sum[(Resnq[i,q] +I*Imsnq[i,q] )*v^(-i),{i,f1}]]

OmegaGQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
f1=Length[Fn];
f=Sum[(ReGnq[i,q] +I*ImGnq[i,q] )(v^(i)+v^(-i)),{i,f1}]]

TFar[q_] := Re[OmegaFar[q]]
TMatrix[q_] := Re[OmegaFar[q]+OmegaSQ[q]]
TQ[q_]:= Re[OmegaGQ[q]]



FnFinal:=Table [If[i < (Length[Fn]+1),Fn[[i]],0.],{i,n}] 
SnqFinal[q_]:=Table [If[i < (Length[Fn]+1),Resnq[i,q] +I*Imsnq[i,q],0.],{i,n}]
GnqFinal[q_]:=Table [If[i < (Length[Fn]+1),ReGnq[i,q] +I*ImGnq[i,q],0.],{i,n}]



(*They are additional terms in the boundary condition of the inclusion caused by temperature effects. They are implemented at Far-field codes.*)

TempratureInteraction1[q_,i_]:=Module[{f,f1,f2,p1,m1,size,fp1,fm1,sp1,sm1,gp1,gm1},
p1= i+1;
m1=i-1;
size=Length[Fn];
If[p1>size,{fp1=0;sp1=0;gp1=0},{fp1=Fn[[i+1]];sp1=SnqFinal[q][[i+1]];gp1=GnqFinal[q][[i+1]]}];
If[m1>size || m1==0,{fm1=0;sm1=0;gm1=0},{fm1=Fn[[i-1]];sm1=SnqFinal[q][[i-1]];gm1=GnqFinal[q][[i-1]]}];
s1=SnqFinal[q][[1]];
(*Print["fp1: ",fp1,"  sp1: ",sp1,"   gp1: ",gp1,"  fm1: ",fm1,"   sm1: ",sm1,"  gm1: ",gm1, "  s1: ", s1];*)
f1=dd[q]/2./i;
f2=Bettal0*((fp1-fm1)+(sm1-sp1)+s1 (-v0[q])^i)+BettalQ[q]*(gm1-gp1);
f=f1*f2]

TempratureInteraction2[q_,i_]:=Module[{f,f1,f2,f3,p1,m1,size,fp1,fm1,sp1,sm1,gp1,gm1},
p1= i+1;
m1=i-1;
size=Length[Fn];
If[p1>size,{fp1=0;sp1=0;gp1=0},{fp1=Fn[[i+1]];sp1=SnqFinal[q][[i+1]];gp1=GnqFinal[q][[i+1]]}];
If[m1>size || m1==0,{fm1=0;sm1=0;gm1=0},{fm1=Fn[[i-1]];sm1=SnqFinal[q][[i-1]];gm1=GnqFinal[q][[i-1]]}];
s1=SnqFinal[q][[1]];
(*Print["fp1: ",fp1,"  sp1: ",sp1,"   gp1: ",gp1,"  fm1: ",fm1,"   sm1: ",sm1,"  gm1: ",gm1, "  s1: ", s1];*)
f1=dd[q]/2./i;
f2=(v0[q]^(2i));
f3=Bettal0*(0(fm1-fp1)*f2+s1 (-v0[q])^i)-BettalQ[q]*(gm1-gp1)*f2;
f=f1*f3]






End[] 
Protect@@Names["Temperature1Inclusion`*"]
EndPackage[] 

BeginPackage["Temperature1Inclusion`",{"geometry`"}] 
EndPackage[] 





































































































