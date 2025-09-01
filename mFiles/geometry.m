(* ::Package:: *)

 BeginPackage[ "geometry`"];
 Unprotect@@Names["geometry`*"];
 ClearAll@@Names["geometry`*"];


RS::usage ="This is for Residual stress calculations. 0 means no RS!";
(*Input-Output- File*)
(*Input-Output-Start*)
n::usage ="number of harmonic terms";
num1::usage ="Crack propagation criterion. Number 1, is the to search on the surface of inclusion for max energy, number 2 uses SIFs.";
num2::usage ="It determines which process has to be used to calculate SIF; 1 is the formula from MM, 2 is the average SIF over short distance.";
ntot::usage ="total number of inclusions in the system";
Nu0::usage ="Nu0 matrix poisson ratio";
Mu0::usage ="Mu0 matrix shear modulus";

(*Far-Field*)
s11::usage ="s11 Farfield stress field";
s12::usage ="s12 Farfield stress field";
s22::usage ="s22 Farfield stress field";

(*Properties of Inclusions*)
Nu::usage ="Poisson ratio for the inclusions (ntot terms, ntot = number of inclusions).";
MuQ::usage ="Shear modulus for the inclusions (ntot terms. Read information about the inclusion from the input file.";
l1::usage ="l1=major axis of elliptical inclusions (ntot terms).";
l2::usage ="l2=Minor axis of elliptical inclusions (ntot terms).";
Theta::usage ="angle of inclination for each inclusion [-90,90]";
Thetaprim::usage ="angle of inclination for each inclusion [0,360]";
zcenter::usage=" (x+Iy): center point of new inclusions(cracks) in the global coordinate system"; 
zcenternewlist1=" center point of new inclusions(cracks) in the global coordinate system"; 

(*Input-Output-End*)

(*Simple Functions - Start*)
dd::usage ="2d is the inter-foci distance of an inclusion";
Dq::usage ="2d is the real part of the inter-foci distance of an inclusion";
dpq::usage ="dpq=dp+dq";
e::usage ="ellipse aspect ratio";
Zeta0::usage ="Zeta0[p]=The matrix-inclusion inerface of the p particle";
v0::usage ="v0=Exp[Zeta0]";
Chi::usage ="Chi=3-4Nu";
MuBar::usage ="relative shear modulus";
Chi0::usage ="Chi0=3-4Nu0 for plain strain and (3-Nu0)/(1+Nu0) for plane stress.";
E0::usage ="Matrix material Young's Modulus in GPa";
EQ::usage ="Inclusions Young's Modulus in GPa";
EndPoint::usage = "EndPoint[p_]= EndPoint of the inclusion"; 
StartPoint::usage = "StartPoint[p_]= StartPoint of the inclusion"; 
lnew1::usage="lnew+delta: lnew1/2 is a vertical distance from the surface of the inclusion at eta0 angle.";

Xi::usage="  = \[Xi] = \[Zeta] + I\[Eta]"; 
XiQ::usage=" \[Xi]q[q] is a local q-th particle in the global z coordinate system. "; 
z::usage=" z= x + iy; coordinate of global system."; 
(*Simple Functions - End*)



(*Read information from the input file.*)
l1temp::usage ="Major axis of elliptical inclusions (ntot terms).";
l2temp::usage ="Minor axis of elliptical inclusions (ntot terms)";
NuQtemp::usage ="Poisson ratio for the inclusions (ntot terms, ntot = number of inclusions )";
MuQtemp::usage ="Shear modulus for the inclusions (ntot terms. ";
tetatemp::usage ="Inclusion inclination angle [-90,90]";
tetatempprim::usage ="Inclusion inclination angle [0,360]";
zcenter1::usage=" zcenter"; 
(*Read information from the input file.- End*)

(*Fracture Energy *)
(*Fracture Energy input - Start *)
G0::usage = "Fracture energy of host material.";
GQ::usage = "Fracture energy of inclusions.";
kQ::usage = "Surface perfectness of inclusions.";
GGB::usage = "Fracture energy of the grain boundaries.";
vertexList::usage = "VGB: Vertices of the grain boundaries.";
BE::usage=" If use finite boundary put 1 in the Fracture energy file, otherwise 0."; 
GB::usage=" If use grain boundary put 1 in the Fracture energy file, otherwise 0."; 
TS::usage=" If use Thermal Stress put 1 in the Fracture energy file, otherwise 0.";
(*ltemp = Read[ss,Number]*) 
lnew::usage="lnew is a half length of new microcrack in the crack propagation."; 
delta::usage="It defines the distance between new added inclusions.";
l2ratio::usage="It is the 1/e ratio, where e is the aspect ratio.";
ttest::usage="Is the number of cracks that are being added to the main crack to predict crack propagation more accuratly.";
inc::usage="";
naccuracy::usage="";

(*Fracture Energy input - End *)
 
GintQ::usage ="Fracture energy of interface of the inclusions (G0*kQ[i]).";
(*Fracture Energy - End *)

(*Crack Propagation- Start*)
lnewtest::usage="This is the obtained function to calculate the optimum distance between two micro-cracks to be able to approximate it with a kinked crack.";
LNEW::usage="LNEW[i]";
LNEW1::usage="LNEW1";
ChangeLNEW::usage="ChangeLNEW[ninitial1_]";
ChangeLNEW1::usage="ChangeLNEW1[deltat_,lt_]";
Changelnewsize::usage=" Changelnewsize[lnewt_,deltat_]";
Deltatest::usage="Deltatest[l2_,l1t_]";
DELTANEW::usage="DELTANEW[initial1_]";
DELTANEW1::usage="DELTANEW[initial1_]";
(*Crack Propagation- End*)

(*Temperature start*)
(*ThermalExpansion input file - Start*)
ThermalExpansion0::usage="ThermalExpansion0= Matrix material thermal expansion 1/K";
ThermalExpansionQ::usage="ThermalExpansionQ= Inclusions' thermzal expansion 1/K";
ThermalConductivity0::usage="Thermal Conductivity of the matrix material. (watts/ meter/ Kelvin)";
ThermalConductivityQ::usage="Thermal Conductivity of the inclusions. (watts/ meter/ Kelvin)";
TStart::usage="Initial temperature! (K)";
TEnd::usage="End temperature! (K) ";
(*ThermalExpansion input file - End*)

ThermalExpansiontemp::usage="An array of Inclusions' thermal expansion 1/K";
ThermalConductivitytemp::usage="Thermal Conductivity of the inclusions. (watts/ meter/ Kelvin)";
ThermalConductivityQBar::usage="Thermal Conductivity of the inclusions/Thermal Conductivity of Matrix material.";
(*Temperature End*)

FractureEnergyBC::usage="FractureEnergyBC[BE1,GB1]";
nChange::usage="nChange[naccuracy_]";

(*Define path of Files and directories*)
FilePath::usage ="Please choose directory of the package.";
OutputFile::usage ="Path to Output file";
InputFile::usage ="Path to input file";
FractureEnergyOutputFile::usage ="Path to Fractureenergy file in output folder";
FractureEnergyInputFile::usage ="Path to Fractureenergy file in input folder";
ThermalExpansionFile::usage ="Path to ThermalExpansion file in output folder";
(*mFilesPath::usage ="Path to mFiles folder.";*)

ReadGeometry::usage ="Read Geometry from the Output file";
FractureEnergy::usage ="Read FractureEnergy file from the Output folder";
ThermalExpansion::usage ="Read ThermalExpansion from the Output file";


(*fg=(Chi0+1)/8/Mu0; (*1/GPa*)
kc=5; (*M N*m^(-3/2)*)
G0=1000*fg*kc^2  (*Pa m = Jul/m^2*)*)


(* ::Input:: *)
(**)


Begin[ "Private`"];
 
(*Define path of Files and directories*)
FilePathPrivet = OpenRead[FileNameJoin[{$HomeDirectory,"DCD_Path"}]];
FilePathPrivet1 = Read[FilePathPrivet];
FilePath := FilePathPrivet1;
(*Print["FilePath = ",FilePath]*)
(*mFilesPath = FileNameJoin[{FilePath, "mFiles"}]*)
InputFile = FileNameJoin[{FilePath, "input", "input"}];
OutputFile = FileNameJoin[{FilePath, "output", "output"}];
FractureEnergyOutputFile = FileNameJoin[{FilePath, "output", "FractureEnergy"}];
FractureEnergyInputFile = FileNameJoin[{FilePath, "input", "FractureEnergy"}];
ThermalExpansionFile  = FileNameJoin[{FilePath, "output", "ThermalExpansion"}];






(*This is for Residual stress calculations. 0 means no RS!*)
RS=0


(*Open and read files!*)


ReadGeometry:= {
s=OpenRead[OutputFile];
n= Read[s,Number];
num1= Read[s,Number];
num2= Read[s,Number];
ntot= Read[s,Number];
Nu0= Read[s,Number];
Mu0= Read[s,Number];
s11= Read[s,Number];
s12= Read[s,Number];
s22= Read[s,Number];
NuQtemp =Flatten[ReadList[s,Expression,1]];
MuQtemp =Flatten[ReadList[s,Expression,1]];
l1temp = Flatten[ReadList[s,Expression,1]];
l2temp = Flatten[ReadList[s,Expression,1]];
tetatemp =Flatten[ReadList[s,Expression,1]];
tetatempprim =Flatten[ReadList[s,Expression,1]];
zcenter1 =Flatten[ReadList[s,Expression,1]];
zcenternewlist1 =Flatten[ReadList[s,Expression,1]];
Close[s];}





FractureEnergy:= { ss=OpenRead[FractureEnergyOutputFile];
G0= Read[ss,Number];
GQtemp = Flatten[ReadList[ss,Expression,1]];
kQtemp = Flatten[ReadList[ss,Expression,1]];
GGBtemp = Flatten[ReadList[ss,Expression,1]];
vertexList = ReadList[ss,Expression,1];
BE = Read[ss,Number];
GB= Read[ss,Number];
TS= Read[ss,Number];
ltemp = Read[ss,Number];
lnew = Read[ss,Number];
delta = Read[ss,Number];
l2ratio = Read[ss,Number];
ttest = Read[ss,Number];
inc = Read[ss,Number];
naccuracy = Read[ss,Number];
Close[ss];
}


ThermalExpansion:= { 
f=OpenRead[ThermalExpansionFile];
ThermalExpansion0 = Read[f,Number];
ThermalExpansiontemp = Flatten[ReadList[f,Expression,1]];
ThermalConductivity0 = Read[f,Number];
ThermalConductivitytemp = Flatten[ReadList[f,Expression,1]];
TStart = Read[f,Number];
TEnd = Read[f,Number];
Close[f];}


ReadGeometry



FractureEnergy



If[TS==1,{
ThermalExpansion;
ThermalExpansionQ[i_]:=Block[{f},f=ThermalExpansiontemp[[i]]];
ThermalConductivityQ[i_]:=Block[{f},f=ThermalConductivitytemp[[i]]];
ThermalConductivityQBar[i_]:=Block[{f,ThermalConductivity0t}, ThermalConductivity0t = If[ThermalConductivity0==0, 0.0000001,ThermalConductivity0];
f = ThermalConductivitytemp[[i]]/ThermalConductivity0t];
}]




Nu[p_] := Block[{f},f=NuQtemp[[p]]]
MuQ[p_] :=Block[{f},f=MuQtemp[[p]]]
(*3 - 4*Nu[q] This is for Plane strain*)
E0 = 2.*Mu0(1+Nu0);
EQ[i_] :=Block[{f},f=2.*MuQ[i](1+Nu[i])] 
Chi[q_] :=Module[{f,f1}, f = (3 - Nu[q])/(1+ Nu[q]);f1= 3 - 4*Nu[q]]
Chi0 = 3 - 4*Nu0;(3 - Nu0)/(1+ Nu0);(*3 - 4*Nu0 This is for Plane strain*)

MuBar[q_] := Module[{f}, f = MuQtemp[[q]]/ Mu0]
l1[p_]:=Block[{f},f=l1temp[[p]]]
l2[p_]:=Block[{f},f=l2temp[[p]]]
e[p_] := Module[{f}, f = l2[p]/l1[p]]
Theta[p_] := Block[{f},f= tetatemp[[p]] Degree]
Thetaprim[p_] := Block[{f},f= tetatempprim[[p]] Degree]
Dq[p_]:=Module[{f},f= Sqrt[l1[p]^2 - l2[p]^2]];
dd[p_] := Module[{f1, f2}, f1 = Sqrt[l1[p]^2 - l2[p]^2];
  f2 = Exp[I Theta[p]] f1]  
dpq[p_, q_] := Module[{f}, f = dd[p] + dd[q]]
Zeta0[p_] := Module[{f}, f = 0.5 Log[(1 + e[p])/(1 - e[p])]]
v0[p_] := Module[{f}, f = Exp[Zeta0[p]]]
Theta[p_,q_]:=Module[{f},f=Theta[q]-Theta[p]]
zcenter[p_] := Block[{f},f= zcenter1[[p]]]
zcenternewtemp[p_] := Block[{f},f= zcenternewlist1[[p]]]
(*v = Exp[Xi]*)
XiQ[q_] := Module[{f, zq}, zq = z - zcenter[q]; f = ArcCosh[zq/dd[q]]]



EndPoint[p_]:=Block[{f},f=zcenter1[[p]]+l1[p]*Exp[I*Thetaprim[p]]]
StartPoint[p_]:=Block[{f},f=zcenter1[[p]]-l1[p]*Exp[I*Thetaprim[p]]]



GQ[p_] := Block[{f},f=GQtemp[[p]]]
kQ[p_] :=Block[{f},f=kQtemp[[p]]]
GintQ[p_] := Block[{f},f=kQtemp[[p]]]
GGB[p_] := Block[{f},f=GGBtemp[[p]]]


(*Define the size of new inclusion in the crack propagation process and the 
distance delta between microcracks.*)

(*This is the obtained function to calculate the optimum distance between 
two micro-cracks to be able to approximate it with a kinked crack.*)

lnewtest[deltat_,l1t_]:=Block[{l2,f,f1,f2,i},f={deltat/l2==0.4917 -(1.455 l2)/l1t+(2.0599 l2^2)/l1t^2-(0.98486 l2^3)/l1t^3};(*{deltat/l2\[Equal]0.531708\[VeryThinSpace]-(2.08999 l2)/l1t+(4.82953 l2^2)/l1t^2-(5.20725 l2^3)/l1t^3+(2.05299 l2^4)/l1t^4};*)
f1={l2}/.Solve[f ,l2];
Print["f1  ",f1];
f2=Flatten[Select[f1,Element[#,Reals]&]];
If[f2=={},lnew,
Table[If[0<f2[[i]]/l1t<=1,f2[[i]],Unevaluated[Sequence[]]],{i,Times@Length[f2]}][[1]]]]

Deltatest[l2_,l1t_]:=Block[{deltat,f,f1,f2,i},f=l2*(0.4917 -(1.455 l2)/l1t+(2.0599 l2^2)/l1t^2-(0.98486 l2^3)/l1t^3)(*{deltat/l2\[Equal]0.531708\[VeryThinSpace]-(2.08999 l2)/l1t+(4.82953 l2^2)/l1t^2-(5.20725 l2^3)/l1t^3+(2.05299 l2^4)/l1t^4};*)]


lnew1=lnew+delta;

LNEW[i_]:=Block[{l,f,lnew0},
Unprotect[lnew,lnew1];
If[i!=1,lnew0=lnewtest[delta,l1[ntot-1]];Clear[lnew,lnew1];lnew=lnew0;
lnew1=lnew+delta,lnew0=lnewtest[delta,l1[ntot]];Clear[lnew,lnew1];lnew=lnew0;
lnew1=lnew+delta]]

LNEW1[initial1_] := Block[{l, f, lnew0, nt}, Unprotect[lnew, lnew1];
  nt = ntot - initial1 - 1;
  lnew0 = lnewtest[delta, l1[ntot - nt]]; Clear[lnew, lnew1]; 
  lnew = 1*lnew0;
  lnew1 = lnew + delta]

LNEW2[deltat_,lt_] := Block[{l, f, lnew0}, Unprotect[lnew, lnew1];
  lnew0 = lnewtest[deltat, lt]; Clear[lnew, lnew1]; 
  lnew = lnew0;
  lnew1 = lnew + delta]

DELTANEW[initial1_] := Block[{l, f,  delta0, nt}, Unprotect[delta, lnew1];
  nt = ntot - initial1 - 1;
 delta0 = Deltatest[lnew, l1[ntot - nt]]; Clear[delta, lnew1]; 
  delta = 1*delta0;
  lnew1 = lnew + delta]

DELTANEW1[l2t_,initial1_] := Block[{l, f,  delta0, nt}, Unprotect[delta, lnew1,lnew];
  nt = ntot - initial1 ;
 delta0 = Deltatest[l2t, l1[ntot - nt]]; Clear[delta, lnew1,lnew]; 
lnew=l2t;
  delta = delta0;
  lnew1 = lnew + delta]


Print["geometry "];
(*LNEW1;*)
ChangeLNEW[ninitial1_,i_] := Block[{s5}, If[i==1,LNEW1[ninitial1],DELTANEW[ninitial1]];
  s5 = OpenWrite[FractureEnergyOutputFile];
  Write[s5, G0];
  Write[s5, GQtemp];
  Write[s5, kQtemp];
  Write[s5, GGBtemp];
  Write[s5, vertexList //. {x_List} :> x];
  Write[s5, BE];
  Write[s5, GB];
  Write[s5, TS];
  Write[s5, ltemp];
  Write[s5, lnew];
  Write[s5, delta];
  Write[s5, l2ratio];
  Write[s5, ttest];
  Write[s5, inc];
  Write[s5, naccuracy];
  Close[s5];
Get["geometry.m"];
Needs["geometry`"]]



ChangeLNEW1[deltat_,lt_] := Block[{s5}, LNEW2[deltat,lt];
  s5 = OpenWrite[FractureEnergyOutputFile];
  Write[s5, G0];
  Write[s5, GQtemp];
  Write[s5, kQtemp];
  Write[s5, GGBtemp];
  Write[s5, vertexList //. {x_List} :> x];
  Write[s5, BE];
  Write[s5, GB];
  Write[s5, TS];
  Write[s5, ltemp];
  Write[s5, lnew];
  Write[s5, delta];
  Write[s5, l2ratio];
  Write[s5, ttest];
  Write[s5, inc];
  Write[s5, naccuracy];
  Close[s5];
Get["geometry.m"];
Needs["geometry`"]]



Changelnewsize[lnewt_,deltat_] := Block[{s5}, 
  s5 = OpenWrite[FractureEnergyOutputFile];
  Write[s5, G0];
  Write[s5, GQtemp];
  Write[s5, kQtemp];
  Write[s5, GGBtemp];
  Write[s5, vertexList //. {x_List} :> x];
  Write[s5, BE];
  Write[s5, GB];
  Write[s5, TS];
  Write[s5, ltemp];
  Write[s5, lnewt];
  Write[s5, deltat];
  Write[s5, l2ratio];
  Write[s5, ttest];
  Write[s5, inc];
  Write[s5, naccuracy];
  Close[s5];
Get["geometry.m"];
Needs["geometry`"]]


FractureEnergyBC[BE1_,GB1_]:={Print["--------------------------------FractureEnergyBC --------------------------------"];
 sFE = OpenWrite[FractureEnergyOutputFile];
  Write[sFE, G0];
  Write[sFE, GQtemp];
  Write[sFE, kQtemp];
  Write[sFE, GGBtemp];
  Write[sFE, vertexList //. {x_List} :> x];
  Write[sFE, BE1];
  Write[sFE, GB1];
  Write[sFE, TS];
  Write[sFE, ltemp];
  Write[sFE, lnew];
  Write[sFE, delta];
  Write[sFE, l2ratio];
  Write[sFE, ttest];
  Write[sFE, inc];
  Write[sFE, naccuracy];
  Close[sFE];
}



(*This function re-write output file with the nacuracy instead of the n number of
 harmunics. It is only to increase n to obtain better SIFs. *)
(*It is not used!*)
nChange[naccuracy_] := {
  s8 = OpenWrite[OutputFile];
  Write[s8, naccuracy];
  Write[s8, num1];
  Write[s8, num2];
  Write[s8, ntot];
  Write[s8, Nu0];
  Write[s8, Mu0];
  Write[s8, s11];
  Write[s8, s12];
  Write[s8, s22];
  Write[s8, NuQtemp];
  Write[s8, MuQtemp];
  Write[s8, l1temp];
  Write[s8, l2temp];
  Write[s8, tetatemp];
  Write[s8, tetatempprim];
  Write[s8, zcenter1];
  Write[s8, zcenternewlist1];
  Close[s8];}


 End[];
 Protect@@Names["geometry`*"]
 EndPackage[]



























