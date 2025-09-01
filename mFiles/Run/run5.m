(* ::Package:: *)

z


BeginPackage["run5`"]
Unprotect@@Names["run5`*"];
ClearAll@@Names["run5`*"];

Pathpoint::usage ="It provides a crack propagation path points.";
l1temp5::usage ="major semi - axes of the ellips read from output1.";
tetatemp5::usage ="Inclusion inclination angle read from output1 file. ";
zcenter5::usage=" center point of new inclusions(cracks) in the global coordinate system from output1 file."; 
plotpath::usage ="Plot crack propagation path.";

Begin["`Private`"] 





Print[Style["There is an issue when the number of initial inclusions is 1, Please check it!",FontSize->20,FontWeight->Bold, FontColor->Red]]
Get["geometry.m"]
(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"];
Needs["geometry`"];



Output1File = FileNameJoin[{FilePath, "output", "output1"}];

s9=OpenRead[Output1File]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
 Read[s9,Number]
NuQtemp5 =Flatten[ReadList[s9,Expression,1]]
MuQtemp5 =Flatten[  ReadList[s9,Expression,1]]
l1temp5 = Flatten[ ReadList[s9,Expression,1]]
l2temp5= Flatten[ ReadList[s9,Expression,1]]
tetatemp5=Flatten[  ReadList[s9,Expression,1]]
tetatemp50=Flatten[  ReadList[s9,Expression,1]]
zcenter5 =Flatten[  ReadList[s9,Expression,1]]
zcenternewlist5 =Flatten[  ReadList[s9,Expression,1]]
Close[s9]

s10=OpenRead[StartFile]
nstart= Read[s10,Number]
 lcrackstart=Read[s10,Number]
Close[s10]

(*nstart = If[nstart0==0, 0, nstart0];*)

Pathpoint:=Block[{r2,r3,r4,r5,r6},
r2=Drop[zcenter5,nstart-1];
r3=Drop[Table[zcenter5[[i]]+(l1temp5 [[i]]) Exp[I*tetatemp5[[i]]Degree],{i,Times@@Dimensions[zcenter5]}],nstart-1];
r4=zcenter5[[nstart]]-(l1temp5 [[nstart]]) Exp[I*tetatemp5[[nstart]]Degree];
r5=Join[r2,Insert[r3,r4,1]];
r6=Sort[Table[{Re[r5[[i]]],Im[r5[[i]]]},{i,Times@@Dimensions[r5]}]]]
(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];

Get["geometry.m"]
Needs["geometry`"]
r5 = Flatten[
  Table[{zcenter1[[p]] -l1[p]*Exp[I*Thetaprim[p]], zcenter1[[p]], 
    zcenter1[[p]] + l1[p]*Exp[I*Thetaprim[p]]}, {p, 
    Times @@ Dimensions[zcenter1]}]]
r6 = Table[{Re[r5[[i]]], Im[r5[[i]]]}, {i, Times @@ Dimensions[r5]}]
ListLinePlot[r6, PlotMarkers -> Automatic,PlotRange\[Rule]Full]*)
plotpath:=Block[{r1,r2},r1=zcenter5 [[nstart]]-(l1temp5 [[nstart]]) Exp[I*tetatemp5[[nstart]]Degree];r2=ListLinePlot[Pathpoint,Frame->True,PlotRange->{Im[r1],6},FrameLabel->{"x(\[Micro]m)","y(\[Micro]m)"}]]              
(*Import["~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/PIC/crackpath.png"]*)

(*
(*This will plot all inclusions in the system!*)
AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];

Get["geometry.m"]
Needs["geometry`"]

rtest[x_,y_]:=Block[{f},f=Table[((x-Re[zcenter[i]])/Cosh[Zeta0[i]])^2 +((y-Im[zcenter[i]])/Sinh[Zeta0[i]])^2 ==dd[i]^2,{i,3}]]
rtest[x,y]

pls=Table[ContourPlot[((x-Re[zcenter[i]])/Cosh[Zeta0[i]])^2 +((y-Im[zcenter[i]])/Sinh[Zeta0[i]])^2 ==dd[i]^2,{x,-10,10},{y,-10,10}],{i,3}]
Show[Table[ContourPlot[((x-Re[zcenter[i]])/Cosh[Zeta0[i]])^2 +((y-Im[zcenter[i]])/Sinh[Zeta0[i]])^2 ==dd[i]^2,{x,-10,10},{y,-10,10}],{i,3}]]

restylePlot2[p_,op:OptionsPattern[ListLinePlot]]:=ListLinePlot[Cases[Normal@p,Line[x__]:>x,\[Infinity]],op,Options[p]]

restylePlot2[pls,Joined\[Rule]True,Axes\[Rule]False,PlotRange\[Rule]{{-10,10},{-10,10}},FrameLabel->{"x(\[Micro]m)","y(\[Micro]m)"},Frame\[Rule]True(*,PlotLegends\[Rule]{"","",""}*)]
*)



End[] 
Protect@@Names["run5`*"]
EndPackage[] 

BeginPackage["run5`",{"geometry`"}] 
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
