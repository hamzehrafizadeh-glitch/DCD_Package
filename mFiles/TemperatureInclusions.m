(* ::Package:: *)

BeginPackage["TemperatureInclusions`"]
Unprotect@@Names["TemperatureInclusionss`*"];
ClearAll@@Names["TemperatureInclusions`*"];



(*In this file thermal stress fields are calculated for circular shape inclusions. 

Then fields in z coordinate system are transformed to the local zq coordinate system 
and then it is rotated to the yq local coordinate system 
(SigmaxyYq, SigmaxxYq, SigmayyYq). 

SEtaEta: is a tangential stress field caused by the grain boundary on the premier of 
the qth inclusion. This function will be used to find a correct propagation direction 
in eta0 function.*)






gamma::usage ="remove it later.";
Bettal0::usage ="For Matrix material: 2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain; 2.*ThermalExpansion0; is for plane stress";
BettalQ::usage ="For Inclusions: 2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain; 2.*ThermalExpansion0; is for plane stress";

Fn::usage ="Fn[n_,q_]: It is a series of the Expansion coefficient of the Far-field complex temperature field.";
Snq::usage ="snq[n_,q_]: The expansion coefficient of the disturbance temperature field induced by the inclusion in the matrix material vanishes at infinity. ";
fnq::usage ="fnq[n_,q_]: It is a series of the Expansion coefficient of the Far-field complex temperature function for many-body problem..";
Gnq::usage ="ReGnq[n_,q_]: The expansion coefficient of the disturbance temperature field in the inclusion. ";
fnqr::usage ="fnqr[n_,q_]: It is a series of the Expansion coefficient of the complex temperature function (caused by other inclusions) for many-body problem..";

OmegaFar::usage ="OmegaFar[q_]: Far-field complex temperature function.";
OmegaSQ::usage ="OmegaSQ[q_]: The disturbance temperature field induced by the inclusion in the matrix material and it is vanishing at infinity.";
OmegaGQ::usage ="OmegaGQ[q_]: The complex temperature function inside the inclusion.";
OmegaFarr::usage ="OmegaFarr[q_]: Far-field complex temperature function without external field.";
OmegaFarrz::usage ="OmegaFarrz[q_]: Far-field complex temperature function without external field and z coordinate system.";
TFar::usage ="Far-field temperature function";
TMatrix::usage ="Temperature field in the matrix material.";
TQ::usage ="Temperature field in the inclusion.";
FTem::usage ="Expansion coefficients of the Far-field temperature. ";
STem::usage ="Expansion coefficients of the disturbance temperature field induced by the inclusion in the matrix material vanishes at infinity. ";

bcTem1nq::usage ="bcTem1nq[n_,q_]: It is an additional term in the boundary condition (1) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";
bcTem1vec::usage ="Vector form of  an additional terms in the boundary condition (1) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";
bcTem2nq::usage ="bcTem2nq[n_,q_]: It is an additional term in the boundary condition (2) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";
bcTem2vec::usage ="Vector form of  an additional terms in the boundary condition (2) of the inclusion caused by temperature effects. They are implemented at Far-field codes.";
TField::usage ="TFiled is the total temperature field caused by the inclusion's thermal interactions. ";
TFieldtotal::usage ="To derive the total temperature field, we must add TField with the external temperature field in the z coordinate system.";


Begin["`Private`"] ;





(*If you need to run the file seperatly, you need to set the path*)
(**AppendTo[$Path,FileNameJoin[{FilePath, "mFiles/Multipole_Method"}]];**)
(*Print["FilePath = ",$Path]*)
Get["geometry.m"];
Needs["geometry`"]
Get["VectorMatrixForms.m"];
Needs["VectorMatrixForms`"]
Get["expansionCoefficients.m"]
Get["abmnpq.m"]
Needs["ExpansionCoefficient`"]
Needs["abmnpq`"]
Print["Temperature_1Inclusion = ",TS];


(*

Clear[VecCreate,BCTem1]
(*VecCreate function is a general function that produces a vector by having func[n,q] and return vector named "funcname" which is a (n*ntot)x1 dimension vector in the below form:
 (\[NoBreak]ReFn1[1,1]
ReFn1[2,1]
ReFn1[3,1]
ReFn1[1,2]
ReFn1[2,2]
ReFn1[3,2]

\[NoBreak]). Here we have 2 inclusions and n-the umber of harmonic terms; n is 3 here.*)
VecCreate[funcname_,func_,n_,ntot_]:= Module[{f,ndim,i,j,ind1},
ndim = n*ntot;
ind1=0;
funcname=ConstantArray[0,{ndim}];
f=
For[i=1,i<ntot+1,i++,
     For[j=1,j<n+1,j ++,
ind1=j+n(i-1);
(*Print["i: ", i," j: ",j, " ind1: ",ind1];*)
funcname[[ind1]]=func[j,i]];
]]

(*Test it
Clear[funcname]
VecCreate[funcname,ReFn1,n,ntot]
MatrixForm[funcname]*)


(*MatrixCreate[matname_,Rmnpq_,n_,ntot_]
This function generates a matrix of (n*ntot)x(n*ntot) dimensions with components of R[m,n,p,q] in the below form.Here ntot=2 and n=2.
(\[NoBreak]Rmnpq[1,1,1,1]	Rmnpq[2,1,1,1]	Rmnpq[3,1,1,1]	Rmnpq[1,1,2,1]	Rmnpq[2,1,2,1]	Rmnpq[3,1,2,1]
Rmnpq[1,2,1,1]	Rmnpq[2,2,1,1]	Rmnpq[3,2,1,1]	Rmnpq[1,2,2,1]	Rmnpq[2,2,2,1]	Rmnpq[3,2,2,1]
Rmnpq[1,3,1,1]	Rmnpq[2,3,1,1]	Rmnpq[3,3,1,1]	Rmnpq[1,3,2,1]	Rmnpq[2,3,2,1]	Rmnpq[3,3,2,1]
Rmnpq[1,1,1,2]	Rmnpq[2,1,1,2]	Rmnpq[3,1,1,2]	Rmnpq[1,1,2,2]	Rmnpq[2,1,2,2]	Rmnpq[3,1,2,2]
Rmnpq[1,2,1,2]	Rmnpq[2,2,1,2]	Rmnpq[3,2,1,2]	Rmnpq[1,2,2,2]	Rmnpq[2,2,2,2]	Rmnpq[3,2,2,2]
Rmnpq[1,3,1,2]	Rmnpq[2,3,1,2]	Rmnpq[3,3,1,2]	Rmnpq[1,3,2,2]	Rmnpq[2,3,2,2]	Rmnpq[3,3,2,2]

\[NoBreak]).
Also,line If[p==q,Break[]] Clear[jnd]; guarantees if both terms are the same,the matrix element becomes zero.*)
Clear[MatrixCreate,matname]
MatrixCreate[matname_,Amnpq_,n_,ntot_]:=Module[{f, ndim,ndim2,i,j,p,q,m,ind,ind2,jnd},
ndim = 1*n*ntot;
     ndim2=ndim/2;
ind=2;
     m=0;
matname=ConstantArray[0,{ndim,ndim}];
For[q=1,q<ntot+1,q++,
     For[j=1,j<n+1,j ++,

For[p=1,p<=ntot,p++,For[i=1,i<n+1,i++,
(*If[p==q,Break[]] Clear[jnd];*)
jnd=n*(p-1)+(i-1)+1;
m+=1;
(*, ind-1,jnd*)
	matname[[ind-1,jnd]]=Amnpq[i,j,p,q]
];];
Clear[ind2];ind2=ind;
Clear[ind];ind=ind2+1;

];];]

(* Test it by
Clear[matname]
MatrixCreate[matname,Amnpq,2,3]
MatrixForm[matname]
Dimensions[matname]*)

(*This is the reverse of the vector function.It allocates vector terms according to the specific array of n*ntot length.*)
VecRev[funcname_,func_,n_,ntot_]:= Module[{f,ndim,i,j,ind1},
ndim = n*ntot;
ind1=0;
f=
For[i=1,i<ntot+1,i++,
     For[j=1,j<n+1,j ++,
ind1=j+n(i-1);
(*Print["i: ", i," j: ",j, " ind1: ",ind1];*)
func[j,i]=funcname[[ind1]]];
]]

(*This is the attempt to solve the temperature field of many body inclusion problems.I HAVE TO WRITE DOWN EQUATIONS IN THE DOCUMENTATION FILE.*)


(*The first step to solve the problem is to write fnq as a function of Snq and from that solve the above equation for Snq.Having Snq,fnq and then Gnq can be calculated.From that,we can use these expansion coefficients of temperature fields and derive the temperature terms of the boundary equations.*)

(**)

[M1Tem1,M2Tem1,BC1, Mtem1,BCTem1,Tem1,Rmnpq,Rmat]

BCTem1[n_,q_]:=-(ThermalConductivityQBar[q]-1)(Fn[n,q]+Conjugate[Fn[n,q]]*v0[q]^(2n));

RSmnpq[m_,i_,p_,q_]:=Block[{f},f= Eta[[m,Abs[i],p,q]]Exp[-I Theta[p,q]]]


M1Tem1[m_,n_,p_,q_]:=Module[{f,f1,f2},
f1=(ThermalConductivityQBar[q]-1)*smnpq[m,n,p,q];
f2=(ThermalConductivityQBar[q]+(v0[q]^(2n) + v0[q]^(-2n))/(v0[q]^(2n) - v0[q]^(-2n)))*KroneckerDelta[m,n]*KroneckerDelta[p,q];
f=f1+f2]

M2Tem1[m_,n_,p_,q_]:=Module[{f,f1,f2},
f1=(ThermalConductivityQBar[q]-1)*v0[q]^(2n)*Conjugate[smnpq[m,n,p,q]];
f2=(2./(v0[q]^(2n) - v0[q]^(-2n)))*KroneckerDelta[m,n]*KroneckerDelta[p,q];
f=f1+f2]*)


Clear[M1Tem1,M2Tem1,BC1, Mtem1,Tem1,Rmnpq,Rmat,bcTem1nq,bcTem2nq,BCTem1]

BCTem1[n_,q_]:=Module[{f},f=-(ThermalConductivityQBar[q]-1)(Fn[n,q]+Conjugate[Fn[n,q]]*v0[q]^(2n))];



M1Tem1[m_,n_,p_,q_]:=Module[{f,f1,f2},
f1=(ThermalConductivityQBar[q]-1)*smnpq[m,n,p,q];
f2=(ThermalConductivityQBar[q]+(v0[q]^(2n) + v0[q]^(-2n))/(v0[q]^(2n) - v0[q]^(-2n)))*KroneckerDelta[m,n]*KroneckerDelta[p,q];
f=f1+f2]

M2Tem1[m_,n_,p_,q_]:=Module[{f,f1,f2},
f1=(ThermalConductivityQBar[q]-1)*v0[q]^(2n)*Conjugate[smnpq[m,n,p,q]];
f2=(2./(v0[q]^(2n) - v0[q]^(-2n)))*KroneckerDelta[m,n]*KroneckerDelta[p,q];
f=f1+f2]


Clear[BCTem,BCTem,M1Tem,M1TemR,M2Tem]
Clear[Ffar]
(*Please note that all the function should be written in funct[n_,q_] form. 
Note the n and q order.*)

(*Please note that the Thermal far-field in the qth coordinate system is: dd[q]*Conjugate[gamma]/2.
This is Conjugate[gamma]zcenter[q] + (v[q]+1/v[q])dd[q]*Conjugate[gamma]/2 in Z coordinate system.

It is important when you plot the field.*)
gamma= 0.*I;30.00+I*0.000;21.213203435596427-21.213203435596427*I;
t0 = 20.; (*Uniform temperature*)
(*This is Fn is zq coordinate system*)
OmegaFar0z[q_] := t0 + Conjugate[gamma]zcenter[q]
Fnz[n_,q_]:=Block[{f},f=KroneckerDelta[n,1]dd[q]*Conjugate[gamma]/2.]

(*This is Fn is yq coordinate system*)
OmegaFar0[q_] := t0 + 0*Conjugate[gamma]zcenter[q];
Fn[n_,q_]:=Block[{f},f=KroneckerDelta[n,1]dd[q]*Conjugate[gamma]/2.]

VecCreate[Ffar,Fn,n,ntot](*NOT SURE WE NEED IT!*)
VecCreate[BCTem,BCTem1,n,ntot]

MatrixCreate[M1Tem,M1Tem1,n,ntot];
MatrixCreate[M2Tem,M2Tem1,n,ntot];
RealMTem = M1Tem+M2Tem;
ImMTem = M1Tem-M2Tem;
RealBCTem = Re[BCTem];
ImBCTem = Im[BCTem];
real = LinearSolve[RealMTem,RealBCTem ]
im = LinearSolve[ImMTem,ImBCTem ]
(*STem is a vector of Snq terms.
To derive Snq,we need a reverse vector operation.*)
STem= real+I*im

Clear[Snq1,Snq,maxn]
maxn=n;
VecRev[STem,Snq1,n,ntot]
Snq[n_,q_]:=If[n>0 && n<= maxn,Snq1[n,q],0]
Clear[Snq1,Snq,maxn]
maxn=n;
VecRev[STem,Snq1,n,ntot]
Snq[n_,q_]:=If[n>0 && n<= maxn,Snq1[n,q],0]


(*Snq[n,q]*)

(*First find fnq, then Gnq eq 2.2.9*)
Clear[RSmnpq,fnq,fnq1]

(*fnq1:=Block[{f,f1,f2,Rmnpq},
MatrixCreate[Rmnpq,smnpq,n,ntot];
f1=Rmnpq.STem;
VecRev[f1,f2,n,ntot];
f=f2]*)

Clear[Rmnpq,f1,fnq1,fTem];
MatrixCreate[Rmnpq,smnpq,n,ntot];
f1=Rmnpq . STem;
VecRev[f1,fnq1,n,ntot];
fTem=f1

(*fnq[n_,q_]:=If[n>0 && n<= maxn,Fn[n,q]+fnq1[n,q],0]*)
(*fnq[n_,q_]:=If[n>0 && n<= maxn,Fn[n,q]+fnq1[n,q],If[n==0,OmegaFar0[q],0]]*)

fnq[n_,q_]:=If[n>0 && n<= maxn,Fn[n,q]+fnq1[n,q],If[n<0 && n>=- maxn,Fn[-n,q]+fnq1[-n,q],If[n==0,OmegaFar0[q],0]]]
fnqr[n_,q_]:=If[n>0 && n<= maxn,fnq1[n,q],If[n<0 && n>=- maxn,fnq1[-n,q],0]]
(*It is not working as f1=Rmnpq.STem; is an array.*)

(*MatrixForm[BC1]
MatrixForm[M1Tem]
MatrixForm[M2Tem]*)

(*Snq[1,1]==STem[[1]]
Snq[2,1]==STem[[2]]
Snq[3,1]==STem[[3]]
Snq[1,2]==STem[[4]]
Snq[2,2]==STem[[5]]
Snq[3,2]==STem[[6]]*)



Gnq[m_,q_]:=(Re[Gnq4[m,q]]+I*Im[Gnq2[m,q]])


Gnq5[m_,q_]:=Module[{f,n},
n=Abs[m];f=(0.5/ThermalConductivityQBar[q])((1+ThermalConductivityQBar[q])fnq[n,q]  +(ThermalConductivityQBar[q]-1)Conjugate[ Snq[n,q]]v0[q]^(-2n)+(ThermalConductivityQBar[q]-1)Conjugate[fnq[n,q]]v0[q]^(-2n))]

Gnq4[m_,q_]:=Module[{f,n},
n=Abs[m];
If[n==0,f= fnq[n,q],f=0.5(2fnq[n,q]+Power[v0[q],-n]*((Conjugate[Snq[n,q]]-Snq[n,q])/((v0[q]^(n)-v0[q]^(-n)))+(Conjugate[Snq[n,q]]+Snq[n,q])/((v0[q]^(n)+v0[q]^(-n)))))];
f]
Gnq3[m_,q_]:=Module[{f,f1,f2,f3,n},
n=Abs[m];
If[n==0,f= Re[fnq[n,q]]+ I*Im[fnq[n,q]]/ThermalConductivityQBar[q],{
f1=(v0[q]^(2n)-v0[q]^(-2n));
f2=fnq[n,q] -Conjugate[Snq[n,q]]/f1 -Snq[n,q]*(v0[q]^(-2n))/f1;
f= f2/ThermalConductivityQBar[q]}];
f]

Gnq2[m_,q_]:=Module[{f,f1,f2,f3,n},
n = Abs[m];
f1=(v0[q]^(2n)-v0[q]^(-2n));
f2=If[f1==0,fnq[n,q],fnq[n,q] -Conjugate[Snq[n,q]]/f1 -Snq[n,q]*(v0[q]^(-2n))/f1];
f= f2/ThermalConductivityQBar[q]]

(*This one is not correct.*)
Gnq1[n_,q_]:=(1+1./ThermalConductivityQBar[q])(fnq[n,q] + Snq[n,q])/2. +(1-1./ThermalConductivityQBar[q])Conjugate[fnq[n,q]]v0[q]^(2n)/2.;


Bettal0 = 2.*ThermalExpansion0;2.*ThermalExpansion0*(1 + Nu0);

BettalQ[q_]:=Module[{f,f1},f1=2.*ThermalExpansionQ[q];f=2.*ThermalExpansionQ[q]*(1 + Nu[q]);f1]

(*2.*ThermalExpansion0*(1 + Nu0) This is for Plane strain
 2.*ThermalExpansion0; is for plane stress *)







(*Complex temperature function*)
(*OmegaFar[q_] Complex far temperature field*)
(*OmegaFar[q_] := Module[{f,f1,i,v},
v = Exp[Xi];

f=OmegaFar0[q]+Sum[fnq[i,q](v^(i)+v^(-i)),{i,n}]]

OmegaSQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];

f=Sum[(Snq[i,q] )*v^(-i),{i,n}]]

OmegaGQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];

f=Sum[(Gnq[i,q] )(v^(i)+v^(-i)),{i,n}]+OmegaFar0 [q]]*)

Unprotect[z,Xi,x,y];
OmegaFar[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
(*This is to eliminate the exponential behaviour of fields.*)
f=Sum[fnq[i,q]v^(i)v^(-(i*i-1)/4),{i,1,n}]+Sum[fnq[i,q]v^(i),{i,-n,0}]]

OmegaSQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];

f=Sum[(Snq[i,q] )*v^(-i),{i,1,n}]]

OmegaGQ[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
(*f=(Gnq[0,q]+Sum[Gnq[i,q]( v^(i)+v^(-i)),{i,n}])*Exp[0*I*Theta[q]]*)
f=Sum[Gnq[i,q]( v^(i)),{i,-n,n}]]

TFar[q_] := Re[OmegaFar[q]]
TMatrix[q_] := Re[(OmegaFar[q]+OmegaSQ[q])]
TQ[q_]:= Re[OmegaGQ[q]]


(*This function calculates all other inclusion's interaction in qth coordinate system. Without far-field.*)
(*I don't think that it is a good idea to use regular pqr fields. The best option is to find 
and sum the singular fields in the z coordinate system.*)
OmegaFarr[q_] := Module[{f,f1,i,v},
v = Exp[Xi];
(*This is to eliminate the exponential behaviour of fields.*)
f=Sum[fnqr[i,q]v^(i)v^(-(i*i-1)/4),{i,1,1}]*Exp[0*I*Theta[q]]+Sum[fnqr[i,q]v^(i),{i,-n,0}]*Exp[0*I*Theta[q]]]

OmegaFarrz[q_] := Module[{f,f1,i,v,Xi},
v = Exp[Xi];
Xi= ArcCosh[(z - zcenter[q])/dd[q]];
f=Sum[fnqr[i,q]v^(i)v^(-(i*i-1)/4),{i,1,1}]*Exp[0*I*Theta[q]]+Sum[fnqr[i,q]v^(i),{i,-n,0}]*Exp[0*I*Theta[q]]]


(*"TFiled" is the total temperature field caused by the inclusion's thermal interactions. 
To derive the total temperature field, we must add TField with the external temperature 
field in the z coordinate system.*)

(*TField1 is the same as the TField function just written in For loop form.*)
TField1:=Module[{f,q,v},
Clear[f];
f=0;
For[q=1,q<=ntot,q++,
Clear[v,Xi,z];
v = Exp[Xi];
  Xi= ArcCosh[(z - zcenter[q])/dd[q]];
 f+=Sum[(Snq[i,q])*v^(-i),{i,1,n}]*Exp[I*Theta[q]]];
f]

TField[n_,ntot_]:=Sum[Sum[(Snq[i,q])*(Exp[ArcCosh[(z - zcenter[q])/dd[q]]])^(-i),{i,1,n}]*Exp[I*Theta[q]],{q,ntot}]
TFieldtotal = TField[n,ntot]+Conjugate[gamma]z +t0


(*They are additional terms in the boundary condition of the inclusion caused by temperature effects. They are implemented at Far-field codes.*)


bcTem1nq[n_,q_]:=Module[{f,f1,f2,f3,f4},
f1=Bettal0*dd[q]/n/2.;
f2=-BettalQ[q]*dd[q]/n/2.;
f3=(fnq[n+1,q]-fnq[n-1,q]) - (Snq[n-1,q]-Snq[n+1,q])- Snq[1,q]((-v0[q])^(n));
f4=Gnq[n+1,q]-Gnq[n-1,q];
f=f1*f3+f2*f4]


bcTem2nq[n_,q_]:=Module[{f,f0,f1,f2,f3,f4},
f0=v0[q]^(2n);
f1=Bettal0*dd[q]/n/2.;
f2=-BettalQ[q]*dd[q]/n/2.;
f3=(fnq[n-1,q]-fnq[n+1,q]) f0 - Snq[1,q]((-v0[q])^(n));
f4=(Gnq[n-1,q]-Gnq[n+1,q]) f0;
f=f1*f3+f2*f4]


VecCreate[bcTem1vec,bcTem1nq,n,ntot]

VecCreate[bcTem2vec,bcTem2nq,n,ntot]



End[] 
Protect@@Names["TemperatureInclusions`*"]
EndPackage[] 

BeginPackage["TemperatureInclusions`",{"geometry`","ExpansionCoefficient`","abmnpq`","VectorMatrixForms`"}] 
EndPackage[] 























































































































