(* ::Package:: *)

BeginPackage["Mergecrack`"]
Unprotect@@Names["Mergecrack`*"];
ClearAll@@Names["Mergecrack`*"];



MergeCrack::usage="MergeCrack[ttest_]: The aim of this part is to merge cracks with an equivalent one.";
MergeCrackPQ::usage="MergeCrackPQ[p_,q_]: The aim of this part is to merge cracks p and q with an equivalent one. You have to use OutputTemp[]after that.";
MergeCrackDouble::usage="";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)

Needs["eta0p`"];


(*The aim of this part is to merge cracks with an equivalent one. 
The "thetaarc" and "centernew" determine the inclination angle and center 
of the new equivalent crack. Function MergeCrack calculates and writes 
the latest system configuration into the temporary file "outputtemp".*)

MergeCrack[ttest_] := {
Get["geometry.m"];
Needs["geometry`"];
zc=zcenter[ntot-ttest+1]+(l1[ntot-ttest+1]) Exp[I*Thetaprim[ntot-ttest+1]];
za=StartPoint[ntot-ttest];
zb=EndPoint[ntot-ttest];
zbt=zcenter[ntot-ttest]+l1[ntot-ttest] Exp[I*Thetaprim[ntot-ttest]];
rca=zc-za;
rc=Abs[rca]/2;
thetacarc=myArcTan[Re[rca],Im[rca]];
thetac=thetacarc*180/Pi;
(*thetac1 =If[thetac>90.0,thetac-180,If[thetac<-90,180+thetac,thetac]];*)
thetac1= ArcTan[Sin[thetacarc]/Cos[thetacarc]]*180.0/Pi;
centernew=zc-rca/2.0;(*=za+rc*Exp[I*thetacarc]=za+rca/2.0*)

tetatemp3=ReplacePart[Delete[tetatemp,{{ntot-ttest}}],ntot-ttest->thetac1];
tetatemp3prim=ReplacePart[Delete[tetatempprim,{{ntot-ttest}}],ntot-ttest->thetac];
l1temp3=ReplacePart[Delete[l1temp,{{ntot-ttest}}],ntot-ttest->rc];
l2temp3=ReplacePart[Delete[l2temp,{{ntot-ttest}}],ntot-ttest->(rc/100)];
NuQ3=ReplacePart[Delete[NuQtemp,{{ntot-ttest}}],ntot-ttest->0];
MuQ3=ReplacePart[Delete[MuQtemp,{{ntot-ttest}}],ntot-ttest->0];
zcentertemp3=ReplacePart[Delete[zcenter1,{{ntot-ttest}}],ntot-ttest->centernew];
zcenternewtemp3=
ReplacePart[Delete[zcenternewlist1,{{ntot-ttest}}],ntot-ttest->centernew];

  OutputtempFile = FileNameJoin[{FilePath, "output", "outputtemp"}];
  
  s5 = OpenWrite[OutputtempFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot - 1];
  Write[s5, Nu0];
  Write[s5, Mu0];
  Write[s5, s11];
  Write[s5, s12];
  Write[s5, s22];
  Write[s5, NuQ3];
  Write[s5, MuQ3];
  Write[s5, l1temp3];
  Write[s5, l2temp3];
  Write[s5, tetatemp3];
 Write[s5, tetatemp3prim];
  Write[s5, zcentertemp3];
  Write[s5, zcenternewtemp3];
  Close[s5];}








MergeCrackPQ[p_,q_] := {
Get["geometry.m"];
Needs["geometry`"];
p1=Max[p,q];q1=Min[p,q];
zc=zcenter[p1]+(l1[p1]) Exp[I*Thetaprim[p1]];
za=StartPoint[q1];
zb=EndPoint[q1];
zbt=zcenter[q1]+l1[q1] Exp[I*Thetaprim[q1]];
rca=zc-za;
rc=Abs[rca]/2;
thetacarc=myArcTan[Re[rca],Im[rca]];
thetac=thetacarc*180/Pi;
(*thetac1 =If[thetac>90.0,thetac-180,If[thetac<-90,180+thetac,thetac]];*)
thetac1= ArcTan[Sin[thetacarc]/Cos[thetacarc]]*180.0/Pi;
centernew=zc-rca/2.0;(*=za+rc*Exp[I*thetacarc]=za+rca/2.0*)

tetatemp3=ReplacePart[Delete[tetatemp,{{q1}}],q1->thetac1];
tetatemp3prim=ReplacePart[Delete[tetatempprim,{{q1}}],q1->thetac];
l1temp3=ReplacePart[Delete[l1temp,{{q1}}],q1->rc];
l2temp3=ReplacePart[Delete[l2temp,{{q1}}],q1->(rc/100)];
NuQ3=ReplacePart[Delete[NuQtemp,{{q1}}],ntot-ttest->0];
MuQ3=ReplacePart[Delete[MuQtemp,{{q1}}],q1->0];
zcentertemp3=ReplacePart[Delete[zcenter1,{{q1}}],q1->centernew];
zcenternewtemp3=
ReplacePart[Delete[zcenternewlist1,{{q1}}],q1->centernew];

  s5 = OpenWrite[OutputtempFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot - 1];
  Write[s5, Nu0];
  Write[s5, Mu0];
  Write[s5, s11];
  Write[s5, s12];
  Write[s5, s22];
  Write[s5, NuQ3];
  Write[s5, MuQ3];
  Write[s5, l1temp3];
  Write[s5, l2temp3];
  Write[s5, tetatemp3];
 Write[s5, tetatemp3prim];
  Write[s5, zcentertemp3];
  Write[s5, zcenternewtemp3];
  Close[s5];}



MergeCrackDouble[ttest_] := {Print["------------------------MergeCrack Double ----------------------------------"];
Get["geometry.m"];
Needs["geometry`"];
zc=zcenter[ntot-ttest+1]+(l1[ntot-ttest+1]) Exp[I*Thetaprim[ntot-ttest+1]];
za=StartPoint[ntot-ttest];
zb=EndPoint[ntot-ttest];
zbt=zcenter[ntot-ttest]+l1[ntot-ttest] Exp[I*Thetaprim[ntot-ttest]];
rca=-zc+zb;
rc=Abs[rca]/2;
thetacarc=myArcTan[Re[rca],Im[rca]];
thetac=thetacarc*180/Pi;
thetac1 =If[thetac>90.0,thetac-180,If[thetac<-90,180+thetac,thetac]];
centernew=Chop[zc+rca/2.0,10^(-4)];(*=za+rc*Exp[I*thetacarc]=za+rca/2.0*)

tetatemp3=ReplacePart[Delete[tetatemp,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->thetac1];
tetatemp3prim=ReplacePart[Delete[tetatempprim,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->thetac];
l1temp3=ReplacePart[Delete[l1temp,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->rc];
l2temp3=ReplacePart[Delete[l2temp,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->(rc/100)];
NuQ3=ReplacePart[Delete[NuQtemp,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->0];
MuQ3=ReplacePart[Delete[MuQtemp,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->0];
zcentertemp3=ReplacePart[Delete[zcenter1,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->centernew];
zcenternewtemp3=
ReplacePart[Delete[zcenternewlist1,{{ntot-ttest},{ntot-ttest-1}}],ntot-ttest-1->centernew];

  s5 = OpenWrite[OutputtempFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot - 2];
  Write[s5, Nu0];
  Write[s5, Mu0];
  Write[s5, s11];
  Write[s5, s12];
  Write[s5, s22];
  Write[s5, NuQ3];
  Write[s5, MuQ3];
  Write[s5, l1temp3];
  Write[s5, l2temp3];
  Write[s5, tetatemp3];
 Write[s5, tetatemp3prim];
  Write[s5, zcentertemp3];
  Write[s5, zcenternewtemp3];
  Close[s5];}



End[] 
Protect@@Names["Mergecrack`*"]
EndPackage[] 

BeginPackage["Mergecrack`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","InterfaceFunc`","EllipseInterfaceCheck`","EnergyTest`","LineInterfaceCheck`","Addone`","CombinMMBE`"}]
EndPackage[] 





























































