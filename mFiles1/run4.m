(* ::Package:: *)

BeginPackage["run4`"]
Unprotect@@Names[" run4`*"];
ClearAll@@Names["run4`*"];
 

(*routput1::usage="routput1[t_,inc_]";*)
(*MergeCrack::usage="MergeCrack[ttest_]";*)
CrackLength::usage="CrackLength";
rGbInt::usage="rGbInt[i_]";
initsub::usage="initsub[] (run4.m)";
doSomeWork1::usage="doSomeWork1[ttest_,inc_,naccuracy_,length,q_,Lnewt_,rgb_]";
export::usage="export[]: Export Crack ****";
exportenergy::usage="exportenergy[]: export Exle file of crackenergy";

Begin["`Private`"] 






(*main-crack-One Directional propagation RHS *)



Print[Style["CHECK FRACTURE ENERGYFILE, IT MAY NOT CONTAIN TS OPTION!",FontSize->20,FontWeight->Bold, FontColor->Red]]

CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/input/input","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];
CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/input/FractureEnergy","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/FractureEnergy"];
CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/input/ThermalExpansion","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/ThermalExpansion"];

AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"]


Get["CombinCDCTEMP0.m"];
Needs["CombinCDCTEMP0`"];


n11=n;
BE1=BE;
GB1=GB;

FractureEnergyBC[0,0];
nChange[1];


Put[ntot-1,l1[ntot],"~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/start"];

Print["Memory in used before"]
Print[MemoryInUse[]];
(*Needs["DisplacementFields`"];*)

Print["Get AddInclusion"]
Get["AddInclusion.m"]
Needs["AddInclusion`"]

Print["Get EnergyTest"]
Get["EnergyTest.m"]
Needs["EnergyTest`"]

Print["Get CombinCDCTEMP"]
Get["CombinCDCTEMP.m"];
Needs["CombinCDCTEMP`"];
(*Get["MergeInclusion.m"]
Needs["MergeInclusion`"]
*)


FractureEnergyBC[BE1,GB1];
nChange[n11];


(*Function r is defined first to copy the output file from first round 
simulation, and then append the coordinate of new inclusion without changing 
their size into the output1 file. Using this function crack propagation path 
can be determined. It is not given an output that the simulation is running 
on, but it will simply determined the new particle's coordinates that are 
added to the initial configuration. 
*)
(*inc is the number of inclusions left from the end of the crack that are not consider before merging.*)
(*routput1[t_,inc_]:={Print["--------------------------routput1 ",i," ----------------------------------"];
If[t==1,CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"],
Print["t is :  ",t];
s3=OpenRead["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];
Read[s3,Number];
Read[s3,Number];
Read[s3,Number];
ntemp=Read[s3,Number];'
Read[s3,Number];
Read[s3,Number];
Read[s3,Number];
Read[s3,Number];
Read[s3,Number];
NuQtemp1=Flatten[ReadList[s3,Expression,1]];
MuQtemp1=Flatten[ReadList[s3,Expression,1]];
l1temp1=Flatten[ReadList[s3,Expression,1]];
l2temp1=Flatten[ReadList[s3,Expression,1]];
tetatemp1=Flatten[ReadList[s3,Expression,1]];
tetatemp1prim=Flatten[ReadList[s3,Expression,1]];
zcenter11=Flatten[ReadList[s3,Expression,1]];
zcenternewlist11=Flatten[ReadList[s3,Expression,1]];
Close[s3];

NuQ2=Join[Drop[NuQtemp1,-1-inc],Take[NuQtemp,-(inc+2)]];
MuQ2=Join[Drop[MuQtemp1,-1-inc],Take[MuQtemp,-(inc+2)]];
l1temp2=Join[Drop[l1temp1,-1-inc],Take[l1temp,-(inc+2)]];
l2temp2=Join[Drop[l2temp1,-1-inc],Take[l2temp,-(inc+2)]];
tetatemp2=Join[Drop[tetatemp1,-1-inc],Take[tetatemp,-(inc+2)]];
tetatemp2prim=Join[Drop[tetatemp1prim,-1-inc],Take[tetatempprim,-(inc+2)]];
zcentertemp2=Join[Drop[zcenter11,-1-inc],Take[zcenter1,-(inc+2)]];
zcenternewtemp2=Join[Drop[zcenter11,-1-inc],Take[zcenter1,-(inc+2)]];
DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];
s4=OpenWrite["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];

Write[s4,n];
Write[s4,num1];
Write[s4,num2];
Write[s4,ntemp+1];
Write[s4,Nu0];
Write[s4,Mu0];
Write[s4,s11];
Write[s4,s12];
Write[s4,s22];
Write[s4,NuQ2];
Write[s4,MuQ2];
Write[s4,l1temp2];
Write[s4,l2temp2];
Write[s4,tetatemp2];
Write[s4,tetatemp2prim];
Write[s4,zcentertemp2];
Write[s4,zcenternewtemp2];
Close[s4];]

}
*)


(*The aim of this part is to merge cracks with an equivalent one. 
The "thetaarc" and "centernew" determine the inclination angle and center 
of the new equivalent crack. Function MergeCrack calculates and writes 
the latest system configuration into the temporary file "outputtemp".*)

(*
MergeCrack[ttest_] := {Print["------------------------MergeCrack ",i," ----------------------------------"];
Get["geometry.m"];
Needs["geometry`"];
zc=zcenter[ntot-ttest+1]+(l1[ntot-ttest+1]) Exp[I*Thetaprim[ntot-ttest+1]];
za=zcenter[ntot-ttest]-l1[ntot-ttest] Exp[I*Thetaprim[ntot-ttest]];
zb=zcenter[ntot-ttest+1]-(delta+l1[ntot-ttest+1]) Exp[I*Thetaprim[ntot-ttest+1]];
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

  s5 = OpenWrite[
    "~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp"];
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

*)




(*Function r2[] reads data from the output file and delete all involved 
inclusions in the crack propagation path except the first two of them. 
It copies output to outputtemp, delete output, and then writes the new 
list to the output file. 
The aim of this part is to recalculate stress fields and SIFs for better 
accuracy using "naccuracy" variable to define the desire number of 
harmonic term.*)

r2[naccuracy_]:={Print["--------------------------------r2 ",i," ----------------------------------"];
s6=OpenRead["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
Read[s6,Number];
NuQtemp4=Flatten[ReadList[s6,Expression,1]];
MuQtemp4=Flatten[ReadList[s6,Expression,1]];
l1temp4=Flatten[ReadList[s6,Expression,1]];
l2temp4=Flatten[ReadList[s6,Expression,1]];
tetatemp4=Flatten[ReadList[s6,Expression,1]];
tetatemp4prim=Flatten[ReadList[s6,Expression,1]];
zcenter14=Flatten[ReadList[s6,Expression,1]];
zcenternewlist14=Flatten[ReadList[s6,Expression,1]];
Close[s6];

r=Table[{ntot-i},{i,0,ntot-ninitial-2}];
nmin=Times@@Dimensions[r];
tetatemp5=Delete[tetatemp4,r];
tetatemp5prim=Delete[tetatemp4prim,r];
l1temp5=Delete[l1temp4,r];
l2temp5=Delete[l2temp4,r];
NuQ5=Delete[NuQtemp4,r];
MuQ5=Delete[MuQtemp4,r];
zcentertemp5=Delete[zcenter14,r];
zcenternewtemp5=Delete[zcenternewlist14,r];

CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp"];
DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];
s7=OpenWrite["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];

Write[s7,naccuracy];
Write[s7,num1];
Write[s7,num2];
(*Write[s7,ntot-ttest+1];*)
Write[s7,ntot-nmin];

Write[s7,Nu0];
Write[s7,Mu0];

Write[s7,s11];
Write[s7,s12];
Write[s7,s22];


Write[s7,NuQ5];
Write[s7,MuQ5];

Write[s7,l1temp5];
Write[s7,l2temp5];
Write[s7,tetatemp5];
Write[s7,tetatemp5prim];
Write[s7,zcentertemp5];
Write[s7,zcenternewtemp5];
Close[s7];
}



(*This function re-load run1.m and SIF.m using output file. Then it 
writes some elements to the cracktemp file to export the crack.csv and 
crack.xlsx files.
The cracktemp contains information about length, angle, SIFs, and 
SIFs-crack for a crack that is approximated with two smaller cracks at 
each simulation step.*)

r3[]:={Print["--------------------------------r3 ",i," ----------------------------------"];
PackageReload["run1`", KillShadowing -> True];
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];
Print["SIFs 2:  ", SIF`k[ntot],"   SIFs 1:  ", SIF`kminus[ntot-1]];
Print["SIFs-crack 2:  ",  SIF`kcrack[ntot,ntot-1],"   SIFs-crack 1:  ", SIF`kcrack[ntot-1,ntot]];

(*Write elements into the crack file*)
PutAppend[List[l1temp[[ntot-1]],l1temp[[ntot]],l2temp[[ntot-1]],l2temp[[ntot]],tetatemp[[ntot-1]],tetatemp[[ntot]],Re[SIF`kminus[ntot-1]],Im[SIF`kminus[ntot-1]],Re[SIF`k[ntot]],Im[SIF`k[ntot]],Re[SIF`kcrackminus[ntot-1,ntot]],Im[SIF`kcrackminus[ntot-1,ntot]],Re[SIF`kcrack[ntot,ntot-1]],Im[SIF`kcrack[ntot,ntot-1]]],"~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/cracktemp"];
}




(*This function calculates the propagated length using the data stored 
in the output1 file.*)

CrackLength:=Block[{s10,f},s10=OpenRead["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
Read[s10,Number];
NuQtemp1=Flatten[ReadList[s10,Expression,1]];
MuQtemp1=Flatten[ReadList[s10,Expression,1]];
l1temp1=Flatten[ReadList[s10,Expression,1]];
l2temp1=Flatten[ReadList[s10,Expression,1]];
tetatemp1=Flatten[ReadList[s10,Expression,1]];
tetatemp1prim=Flatten[ReadList[s10,Expression,1]];
zcenter11=Flatten[ReadList[s10,Expression,1]];
zcenternewlist11=Flatten[ReadList[s10,Expression,1]];
Close[s10];
ninitial1=If[ninitial==0,1,ninitial];
f=2.0*Sum[l1temp1[[i]],{i,ninitial1,Times@@Dimensions[l1temp1]}]]



(*This deletes output, copy outputtemp to output, and then deletes outputtemp.*)
(*r4[]:={Print["--------------------------------r4 ",i," ----------------------------------"];
DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];
ClearSystemCache[];
CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output"];
DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp"];};
*)





(*It exports the cracktemp text file to the crack.csv and crack.xlsx files.*)

export[]:={s8 = OpenRead[
   "~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/cracktemp"];
rtest = Table[Flatten[ReadList[s8, Expression, 1]], {j, i+1}];
Close[s8];
Clear[L11, L12, L21, L22, theta1, theta2, kminus1I, k2I, kcrack1I, kcrack2I,kminus1II, k2II, kcrack1II, kcrack2II];
rtest1 = Join[{{L11, L12, L21, L22, theta1, theta2, kminus1I,kminus1II, k2I,k2II, kcrack1I,kcrack1II, kcrack2I,kcrack2II}},rtest];

Export["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/crack.csv",rtest1];
Export["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/crack.xlsx",rtest1];
}


exportenergy[] := {ss8 = 
   OpenRead[
    "~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/crackenergy"];
  rtest2 = Table[Flatten[ReadList[ss8, Expression, 1]], {j, i + 1}];
  Close[ss8];
  rtest3 = Join[{{kI, kII, gc, g, Crack - Length}}, rtest2];
  Export["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/Energy.xlsx", rtest3];}



rGbInt[i_]:={Clear[rGB1,rGB2 ,rInt1,rInt2];
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
rInt=Flatten[Testinterface1[ntot]];
rInt1=rInt[[1]];
rInt2=If[Length[rInt]!=1,rInt[[2]],0];

Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];

rGB=Flatten[LineIntersectCheck[i,ntot]];
rGB1=rGB[[1]];
rGB2=If[Length[rGB]!=1,rGB[[2]],0];

Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
penetration=IntersectCheckQ[ntot];
If[Length[penetration]!=1,rGB1=0;rGB2=0;rt=Flatten[TestPenetration[ntot]];
If[rt[[1]]==2&&rt[[2]]==-1,Print["Penetration -1 & I have to add inclusion"];]];

Print["rGB   ","rGB1   ","rGB2    ","rInt   ","rInt1    ","rInt2   "];
Print[rGB,"     ",rGB1,"      ",rGB2 ,"      ",rInt,"       ",rInt1,"      ",rInt2];
Put[rGB1,rGB2 ,rInt1,rInt2,"~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/crack"];
Clear[rGB1,rGB2 ,rInt1,rInt2];}


(*This is the main function in the run.m file. 
The function recalculates the new configuration of the system  while the 
crack length reaches a desired length. It can eventually predict crack 
propagation path.
 The "ncrack" parameter is a number that determines the number of 
involved inclusions in the propagated crack.*)

doSomeWork1[ttest_,inc_,naccuracy_,length_,q_,Lnewt_,rgb_]:={If[inc>= ttest-1, Print["Oops inc is larger than ttest-2!"];Quit[]];

i=0;While[ttest!=0,i++;Print["doSomework1 i: ",i];

If[i==1,t=ttest,t=1];

Print["-------------CHECK FIRST IF THE CRACK PROPAGATE? -------------------------"];

If[i==1,If[num2==1,PropagationTest[ntot,i,ttest]];If[rgb==1,rGbInt[1]]];
Print["-------------------AddInc-MergeInc ",i," ----------------------------------"];

Get["geometry.m"];
Needs["geometry`"];
If[i!=1,ChangeLNEW[ninitial,Lnewt]];
Print["lnew =  ", lnew, "   delta =  ",delta];
PutAppend[List[lnew,delta],"~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/lnew"];

Get["eta0p.m"];
AddInclusion[t,q];

Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
(*PropagationTest[ntot,2,ttest];*)


(*MergeInclusion[ntot-inc,ntot-inc-1];*)

(*Write output1*)

(*If[tt==1,CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"],r[t]];*)
Print["--------------------------routput1 i ----------------------------------"];
If[i!=1 && Abs[tetatempprim[[ntot]]- tetatempprim[[ntot-1]]]>=90 , Print[Style["Delta theta is Larger than 90 Degree, it cause an uncertainty in the results!", FontColor -> Purple]]];
routput1[i,inc,1];

(*Change the output to the new output with a crack. The crack is consistes of two smaller crack*)
(*r2[naccuracy];

r3[];
(*delet output, copy outputtemp to output, then delet outputtemp*)
r4[];*)

Print["------------------------MergeCrack ",i," ----------------------------------"];
MergeCrack[ttest];
Print["CrackLength/lltest:  ", CrackLength/lltest];
Print["Crack Length is:  ",CrackLength,"; lltest:  ",lltest];

Print["--------------------------------r4 ",i," ----------------------------------"];
OutputTemp[];


PassInclusion[ntot-ttest];

If[(CrackLength/lltest)> length,(*Get["geometry.m"];
Needs["geometry`"];r2[naccuracy];r3[];(*For[kt=ninitial,kt<=IntegerPart[ntot/2]+1,kt++,MergeFinal[kt+1,kt]];
Get["geometry.m"];
Needs["geometry`"];
r3[];
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];*)

*)
Break[]];
]
(*The file contains below information.
l1[1],l1[2],l2[1],l2[2],kminus[1],theta[1],theta[2],k[2],kcrack[1],kcrack[2]*)
(*export[];*)
exportenergy[];

(*This file produces crack path!*)
Get["run5.m"];
Needs[ "run5`"];
}







Needs["SubKernels`LocalKernels`"]
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]])
initsub[]

(*Total file runs from here! *)
(*Evaluate[doSomeWork1[4,0,3],subkernel]/.{$Failed:>(Print@"Ouch!";initsub[])}
Print[MemoryInUse[]];
Print["Time in used"]
Print[TimeUsed[]];*)



End[] 
Protect@@Names["run4`*"]
EndPackage[] 

(*BeginPackage["run4`",{ "run5`" }] 
EndPackage[] *)
(*, "run5`"*)
BeginPackage["run4`",{"run1`","PackageManipulations`","EnergyTest`"}] 
EndPackage[] 













(* ::Input:: *)
(**)
