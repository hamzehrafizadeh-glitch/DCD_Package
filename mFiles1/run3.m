(* ::Package:: *)

BeginPackage["run3`"]
Unprotect@@Names["run3`*"];
ClearAll@@Names["run3`*"];
Begin["`Private`"] 






(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"]
Needs["geometry`"]
Get["PackageManipulations.m"]
Needs["PackageManipulations`"]

DeleteFile[FileNameJoin[{FilePath, "output", "outputfinal"}]]
CopyFile[InputFile,OutputFile]
CopyFile[InputFile,FileNameJoin[{FilePath, "input", "input-initial"}]]

Print["Memory in used before"]
Print[MemoryInUse[]];

t = 7;
Clear[j]
j=0;
doSomeWork[t_] := {For[q = 1, q < t + 1, q++, Print["doSomework q: ",q];
   If[q == 1,
    Get["run1.m"];
    Needs["run1`"];
    Get["run2.m"]; Needs["run2`"]; Print["ii ", run2`ii]; 
    If[run2`ii == 1, Break[]];
	j++;
	Print["j is:  ", j],
    PackageReload["run1`", KillShadowing -> True];
    CloseKernels[subkernel]; initsub[];
    PackageReload["run2`", KillShadowing -> True];
    (*CloseKernels[subkernel];initsub[]*)Print["ii: ", run2`ii];
    If[run2`ii == 1, Break[]];
	j ++;
    Print["j is:  ", j];
    ]]}



Needs["SubKernels`LocalKernels`"]
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]])
initsub[]

Evaluate[doSomeWork[t],subkernel]/.{$Failed:>(Print@"Ouch!";initsub[])}
Print[MemoryInUse[]];
Print["Time in used"]
Print[TimeUsed[]];


Needs["geometry`"]





CopyFile[OutputFile,FileNameJoin[{FilePath, "output", "outputfinal"}]];
DeleteFile[OutputFile];
(*DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp"];*)



End[] 
Protect@@Names["run3`*"]
EndPackage[] 

BeginPackage["run3`",{"run1`", "run2`","PackageManipulations`" }] 
EndPackage[] 















(* ::Input:: *)
(*(*(*The aim of this part is to merge cracks with an equivalent one. The "thetaarc" and*)
(* "centernew" determine the inclination angle and center of the new equivalent crack.*)*)
(**)
(*zc=zcenter[ntot-ttest+1]+(l1[ntot-ttest+1]) Exp[I*Theta[ntot-ttest+1]]*)
(*za=zcenter[ntot-ttest]-l1[ntot-ttest] Exp[I*Theta[ntot-ttest]]*)
(*zb=zcenter[ntot-ttest+1]-(delta+l1[ntot-ttest+1]) Exp[I*Theta[ntot-ttest+1]];*)
(*zbt=zcenter[ntot-ttest]+l1[ntot-ttest] Exp[I*Theta[ntot-ttest]];*)
(*rca=zc-za*)
(*rc=Abs[rca]/2*)
(*thetacarc=ArcTan[Im[rca]/Re[rca]]*)
(*thetac=ArcTan[Im[rca]/Re[rca]]*180/Pi*)
(*centernew=zc-rca/2.0(*=za+rc*Exp[I*thetacarc]=za+rca/2.0*)*)
(**)
(*tt=1;*)
(*(*Function r is defined first to copy the output file from first round *)
(*simulation, and then append the coordinate of new inclusion without changing *)
(*their size into the output1 file. Using this function crack propagation path *)
(*can be determined. It is not given an output that the simulation is running *)
(*on, but it will simply determined the new particle's coordinates that are *)
(*added to the initial configuration. *)
(**)*)
(*r[]:={s3=OpenRead["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*Read[s3,Number];*)
(*NuQtemp1=Flatten[ReadList[s3,Expression,1]];*)
(*MuQtemp1=Flatten[ReadList[s3,Expression,1]];*)
(*l1temp1=Flatten[ReadList[s3,Expression,1]];*)
(*l2temp1=Flatten[ReadList[s3,Expression,1]];*)
(*tetatemp1=Flatten[ReadList[s3,Expression,1]];*)
(*zcenter11=Flatten[ReadList[s3,Expression,1]];*)
(*zcenternewlist11=Flatten[ReadList[s3,Expression,1]];*)
(*Close[s3];*)
(**)
(**)
(*tetatemp2=Join[tetatemp1,Table[tetatemp[[i]],{i,ntot-t+1,ntot}]];*)
(*l1temp2=Join[l1temp1,Table[l1temp[[i]],{i,ntot-t+1,ntot}]];*)
(*l2temp2=Join[l2temp1,Table[l2temp[[i]],{i,ntot-t+1,ntot}]];*)
(*NuQ2=Join[NuQtemp1,Table[NuQtemp[[i]],{i,ntot-t+1,ntot}]];*)
(*MuQ2=Join[MuQtemp1,Table[MuQtemp[[i]],{i,ntot-t+1,ntot}]];*)
(*zcentertemp2=Join[zcenter11,Table[zcenter11[[i]],{i,ntot-t+1,ntot}]];*)
(*zcenternewtemp2=Join[zcenternewlist11,Table[zcenternewlist1[[i]],{i,ntot-t+1,ntot}]];*)
(**)
(*DeleteFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];*)
(*s4=OpenWrite["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"];*)
(**)
(*Write[s4,n];*)
(*Write[s4,n1];*)
(*Write[s4,n2];*)
(*Write[s4,ntot];*)
(**)
(*Write[s4,Nu0];*)
(*Write[s4,Mu0];*)
(**)
(*Write[s4,s11];*)
(*Write[s4,s12];*)
(*Write[s4,s22];*)
(**)
(**)
(*Write[s4,NuQ2];*)
(*Write[s4,MuQ2];*)
(**)
(*Write[s4,l1temp2];*)
(*Write[s4,l2temp2];*)
(*Write[s4,tetatemp2];*)
(*Write[s4,zcentertemp2];*)
(*Write[s4,zcenternewtemp2];*)
(*Close[s4];*)
(*}*)
(*If[tt==1,CopyFile["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output","~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/output1"],r[]]*)
(**)
(*ttest=3;*)
(*tetatemp3=ReplacePart[Delete[tetatemp,{{ntot-ttest}}],ntot-ttest->thetac]*)
(*l1temp3=ReplacePart[Delete[l1temp,{{ntot-ttest}}],ntot-ttest->rc]*)
(*l2temp3=ReplacePart[Delete[l2temp,{{ntot-ttest}}],ntot-ttest->(rc/100)]*)
(*NuQ3=ReplacePart[Delete[NuQtemp,{{ntot-ttest}}],ntot-ttest->0]*)
(*MuQ3=ReplacePart[Delete[MuQtemp,{{ntot-ttest}}],ntot-ttest->0]*)
(*zcentertemp3=ReplacePart[Delete[zcenter1,{{ntot-ttest}}],ntot-ttest->centernew]*)
(*zcenternewtemp3=*)
(*ReplacePart[Delete[zcenternewlist1,{{ntot-ttest}}],ntot-ttest->centernew]*)
(**)
(*s5=OpenWrite["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/output/outputtemp"];*)
(**)
(*Write[s5,n];*)
(*Write[s5,n1];*)
(*Write[s5,n2];*)
(*Write[s5,ntot-1];*)
(**)
(*Write[s5,Nu0];*)
(*Write[s5,Mu0];*)
(**)
(*Write[s5,s11];*)
(*Write[s5,s12];*)
(*Write[s5,s22];*)
(*Write[s5,NuQ3];*)
(*Write[s5,MuQ3];*)
(**)
(*Write[s5,l1temp3];*)
(*Write[s5,l2temp3];*)
(*Write[s5,tetatemp3];*)
(*Write[s5,zcentertemp3];*)
(*Write[s5,zcenternewtemp3];*)
(*Close[s5];*)*)
