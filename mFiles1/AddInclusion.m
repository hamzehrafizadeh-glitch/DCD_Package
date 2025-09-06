(* ::Package:: *)

BeginPackage["AddInclusion`"]
Unprotect@@Names["AddInclusion`*"];
ClearAll@@Names["AddInclusion`*"];



(*main-crack-One Directional propagation RHS 

AddInclusion[t]
The "AddInclusion" function adds "t" number of inclusions with a specific length "lnew1" to the system. Adding process is started from the last inclusion "ntot-th" in the input file. 
*)


AddOne1::usage="AddOne1[]: If just inclusions interface";
AddOne2::usage="AddOne2[]: If just Grain boundary";
AddOne3::usage="AddOne3[]: If both inclusions interface and GB";
Interface1::usage="Interface1[i_,j_] : i is the intersected inclusion.";
Interface2::usage="Interface2[i_,j_] : i is the grain boundary.";
Interface3::usage="Interface3[i_,j_,p_] : i is the grain boundary, j is the intersected inclusion, and p is the number of last crack ntot.";
AddInclusion::usage="AddInclusion[n_,i_]: i=1) No GB but Interface, i=2) GB, i=3) GB & Interface";
AddInclusionDouble::usage="AddInclusionDouble[n_]";
gc::usage = "Maximum Fracture energy at the crack tipe depending on where the crack is. ";



Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
(*Get["DisplacementFields.m"];
Needs["DisplacementFields`"];*)

(*Get["linearSolve.m"];*)
Print["AddInclusion"]
Needs["eta0p`"];
Needs["SIF`"];
If[BE==1,Needs["run02`"]];
Needs["EllipseInterfaceCheck`"];
Needs["EnergyTest`"];
Needs["LineInterfaceCheck`"];
Needs["Addone`"];
Needs["InterfaceFunc`"];
Needs["Changesize`"];





InterfaceFractureCoeff[p_,ttest_] := Block[{r1, r, r2,f,f1,f2,j,jj},
f  = If[p>ttest,f=Table[j,{j,p-ttest,p-1,1}];
Print["coefficient test=  ", f];
f1 = Sum[means22[p, jj]*Sqrt[Pi*l1[p]],{jj,f}];
f2 = Sum[means12[p, jj]*Sqrt[Pi*l1[p]],{jj,f}];
r = f1 + I*f2 + k[p];
r1 = ArcTan[Im[r]/Re[r]];
r2 = 1 + 0.(Im[r]/Re[r])^2,1]]


(*Needs["SubKernels`LocalKernels`"]
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]])
initsub[]*)




Interface1[i_,j_]:=Block[{GammaInt,GammaFrac,GammaInc,rt,r2,jj,j1},Unprotect[gc];
gc = If[i != 0,GintQ[i], G0];
Print["gc   ",gc];
GammaInt= GintQ[i];
GammaFrac=  G0;
GammaInc=GQ[i];
j1=j;
If[j1==2,
Reload;Print["Reload"],
If[j1==0,
AddOne[ntot];
Get["geometry.m"];
Needs["geometry`"];
rt=Flatten[TestPenetration[ntot]];
Print["Penetration: ", rt];
If[rt[[1]]==2&&rt[[2]]==-1,{Reload;Clear[j1];j1=2;},Unevaluated[Sequence[]]];
Print["r1 = ",j1],r2=AddCrackInterface1[i,GammaInt,GammaFrac,GammaInc];Print["AddInterfacecrack"]]];
jj=If[j1!=2,3,2]]


Interface2[i_,j_]:=Block[{GammaGB,GammaFrac,rt,r2,jj,j1,ii},Unprotect[gc];
ii=1;
gc = G0;
Print["gc   ",gc];
GammaGB= GGB[i];
GammaFrac=  G0;
Print["GammaGB =  ",GammaGB,"  GammaFrac =  ",GammaFrac];
j1=j;
If[j1==2,
Reload;Print["Reload"],
If[j1==0,
AddOne[ntot];
Get["geometry.m"];
Needs["geometry`"];
rt=Flatten[TestPenetration[ntot]];
Print["Penetration: ", rt];
If[rt[[1]]==2&&rt[[2]]==-1,Reload;Clear[j1];j1=2,Unevaluated[Sequence[]]];
Print["j = ",j],
r2=AddGBCrack[i,ntot,GammaGB,GammaFrac];CrackPenetrateGB[GammaFrac,ii];Print["AddGBCrack"]]];
jj=If[j1!=2,3,2]]

(*-----------------------------------Interface3 ----- --------------------------*)
Interface3[i_,j_,p_]:=Block[{jj,f,ThetaNew,t,GammaGB,GammaFrac,GammaInt,rgb1,rgb2,rgb3,t2,rinc1,rinc2,rinc3},
Print["------------------------------Interface3 GB --------------------------"];
Get["geometry.m"];
Needs["geometry`"];
f=GBLine[i];
ThetaNew=f[[5]]*180/Pi;
t=EndPoint[p];
GammaGB= GGB[i];
GammaFrac=  G0;
GammaInt= GintQ[j];
Print["GammaGB =  ",GammaGB,"  GammaFrac =  ",GammaFrac, "  GammaInt =  ", GammaInt];

     rgb1=NewInc[p,ThetaNew,GammaGB,GammaFrac];
     rgb2 =rgb1[[1]];
     rgb3=If[Length[rgb1]!=1,rgb1[[3]]];

Print["---------------------------Interface3 Ellipse -----------------------"];
t2=yprime[j,Re[t],Im[t]];
rinc1=NewInc[p,t2,GammaInt,GammaFrac];
rinc2=rinc1[[2]];
rinc3 =rinc1[[3]];

If[rgb3>rinc3,jj=Interface2[i,1](*;Get["geometry.m"];
Needs["geometry`"];
rt=Flatten[TestPenetration[ntot]];
Print["Penetration: ", rt];
If[rt[[1]]==2&&rt[[2]]==-1,jj=2]*),jj=Interface1[j,1]];
jj]


AddOne1[]:=Block[{r,r1,r2,rt,jj,i,ii,GammaInt,GammaFrac,GammaInc,rr2,rr3},
Reload;
PropagationTest[ntot,2,ttest];

jj=2;While[jj==2,
If[ninitial==0,AddOne[ntot];Clear[jj];jj=3,
r=Flatten[Testinterface1[ntot]];
r1=r[[1]];
i=If[Length[r]!=1,r[[2]],0];
(*Print[Style["ChangeLNEW2  ", FontColor -> Blue], ntot];ChangeLNEW2[i,ntot-ttest-1,delta];*)
jj=Interface1[i,r1]]]]

AddOne2[]:=Block[{r,r1,r2,jj,i,GammaGB,GammaFrac},
(*In this example i=1, just one Grain boundary.*)
i=1;
Reload;
PropagationTest[ntot,2,ttest];

jj=2;While[jj==2,
r=Flatten[LineIntersectCheck[i,ntot]];
r1=r[[1]];
r2=If[Length[r]!=1,r[[2]],0];
jj=Interface2[i,r1]]]




AddOne3[]:=Block[{r,r1,r2,rt,jj,iGB,GammaGB,GammaFrac},
(*In this example i=1, just one Grain boundary.*)
iGB=1;

jj=2;While[jj==2,

s6=OpenRead[FileNameJoin[{$Path,"Crack"}]];
rGB1=Read[s6,Number];
rGB2=Read[s6,Number];
rInt1=Read[s6,Number];
rInt2=Read[s6,Number];
Close[s6];

Reload;
PropagationTest[ntot,2,ttest];
(*If[rInt2!=0,Print[Style["ChangeLNEW2  ", FontColor -> Blue],ntot, rInt2];ChangeLNEW2[rInt2,ntot-ttest-1,delta];Get["geometry.m"],Unevaluated[Sequence[]]];*)
(*CHECK THIS ONE SEEMS RONG TO ME! rInt1==2*)

If[rInt1==2&&rGB1==0,jj=Interface1[rInt2,rInt1]];
If[rGB1==2 && rInt1==0,jj=Interface2[iGB,rGB1]];
If[rInt1==0 && rGB1==0,jj=Interface1[rInt2,rInt1]];
If[rInt1==1 && rGB1==0,jj=Interface1[rInt2,rInt1]];
If[rInt1==0 && rGB1==1,jj=Interface2[iGB,rGB1]];
If[rInt1==1&&rGB1==1,jj=Interface3[iGB,rInt2,ntot]];
If[(rInt1==1&&rGB1==2)||(rInt1==2&&rGB1==2),jj=Interface3[iGB,rInt2,ntot]];
If[rInt1==2&&rGB1==1, jj=2];

Print["------------------------- GB & Int ------------------------------"];
Clear[rGB1,rGB2 ,rInt1,rInt2];
Get["geometry.m"];
Needs["geometry`"];
rInt=Flatten[Testinterface1[ntot]];
rInt1=If[rInt[[1]]!=0,1,rInt[[1]]];
rInt2=If[Length[rInt]!=1,rInt[[2]],0];

Get["geometry.m"];
Needs["geometry`"];

rGB=Flatten[LineIntersectCheck[iGB,ntot]];
rGB1=rGB[[1]];
rGB2=If[Length[rGB]!=1,rGB[[2]],0];
If[rGB1==2 && rGB2==-1, jj=2];

Get["geometry.m"];
Needs["geometry`"];
penetration=IntersectCheckQ[ntot];
Print["penetration =   ",penetration];
If[Length[penetration]!=1,rGB1=0;rGB2=0;rt=Flatten[TestPenetration[ntot]];
If[rt[[1]]==2&&rt[[2]]==-1,Print["Penetration -1 & I have to add inclusion"];]];

Print["rGB   ","rGB1   ","rGB2    ","rInt   ","rInt1    ","rInt2   "];
Print[rGB,"     ",rGB1,"      ",rGB2 ,"      ",rInt,"       ",rInt1,"      ",rInt2];
Put[rGB1,rGB2 ,rInt1,rInt2,CrackOutputFile];
Clear[rGB1,rGB2 ,rInt1,rInt2];
Print["jj =   ", jj];
]]


(*AddInclusion[t_]:={For[q=1,q<t+1,q++,Print["AddInclusion number: ",q];
AddOne3[];
;]}*)


AddInclusion[t_,i_]:={For[q=1,q<t+1,q++,Print["AddInclusion number: ",q];
<<Utilities`CleanSlate`;
CleanSlate[];
ClearInOut[];
(*If just inclusions interface, use AddOne1[];*)
(*If just Grain boundary, use AddOne2[];*)
(*If both inclusions interface and GB, use AddOne3[];*)
If[i==1,AddOne1[]];
If[i==2,AddOne2[]];
If[i==3,AddOne3[]];

]}



AddInclusionDouble[t_]:={For[q=1,q<t+1,q++,Print["AddInclusion number: ",q];
Reload;
AddOneDouble[ntot];
;]}


End[] 
Protect@@Names["AddInclusion`*"]
EndPackage[] 

BeginPackage["AddInclusion`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","InterfaceFunc`","EllipseInterfaceCheck`","EnergyTest`","LineInterfaceCheck`","Addone`","CombinMMBE`"}] 
EndPackage[] 





















































































