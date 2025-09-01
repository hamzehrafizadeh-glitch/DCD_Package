(* ::Package:: *)

BeginPackage["Changesize`"]
Unprotect@@Names["Changesize`*"];
ClearAll@@Names["Changesize`*"];




ChangeSize::usage="ChangeSize[i_, lnew_, centernew_,sign_,tetatnew_,tetatnewtemp_] ";
PlaceiToj::usage="PlaceiToj[i_,j_]; delete inclusion i and put it at j place.";
Replace1::usage="Replace1[i_]: Replace the ith inclusion as the ntot-th crack.";
Replaceij::usage="Replaceij[i_,j_]: Replace the ith inclusion with jth one.";
ReplaceInc::usage="ReplaceInc[j_,ttest_]: Replace the ith inclusion to be part of the crack.";
DeleteInc::usage="DeleteInc[i_]";
swap::usage="swap[a_,i_,j_]: it awaps elemt i and jth of list a.";
OutputTemp::usage="OutputTemp[]: This deletes output, copy outputtemp to output, and then deletes outputtemp.";
Copyoutput::usage="Output: This copy output to outputtemp.";
routput1::usage="routput1[t_,inc_]";
Length2Crack::usage="Length2Crack[naccuracy_]: it reads data from the output file and delete all involved 
inclusions in the crack propagation path except the first two of them. 
It copies output to outputtemp, delete output, and then writes the new 
list to the output file. 
The aim of this part is to recalculate stress fields and SIFs for better 
accuracy using naccuracy variable to define the desire number of 
harmonic term.";
ExportSIF2::usage="ExportSIF2[p_,q_]: This function re-load run1.m and SIF.m using output file. Then it 
writes some elements to the cracktemp file to export the crack.csv and 
crack.xlsx files.
The cracktemp contains information about length, angle, SIFs, and 
SIFs-crack for a crack that is approximated with two smaller cracks at 
each simulation step.";




Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)

Needs["eta0p`"];


ChangeSize[i_, lnew_, centernew_,sign_,tetatnew_,tetatnewtemp_] := 
 Block[{l1temp4, l2temp4, zcentertemp3, s5,tetatemp4,tetatempprim4 ,MuQ1,NuQ1}, 
  Print["----------------Change the ", i, "  th Crack Size----------------"];
  (*Get["geometry.m"];
  Needs["geometry`"];*)
If[sign==0,
{NuQ1=NuQtemp;
MuQ1=MuQtemp;
l1temp4 = ReplacePart[l1temp, i -> lnew];
l2temp4 = ReplacePart[l2temp, i ->lnew/l2ratio];
tetatemp4 = ReplacePart[tetatemp, i ->tetatnew];
tetatempprim4 = ReplacePart[tetatempprim, i ->  tetatnewtemp];
zcentertemp3 = ReplacePart[zcenter1, i -> centernew];},
If[sign==-1,
{NuQ1=Delete[NuQtemp,i];
MuQ1=Delete[MuQtemp,i];
l1temp4 = ReplacePart[Delete[l1temp, i ], i-1 -> lnew];
l2temp4 = ReplacePart[Delete[l2temp, i ], i-1 ->lnew/l2ratio];
tetatemp4 =ReplacePart[Delete[tetatemp, i ], i-1 ->tetatnew];
tetatempprim4 = ReplacePart[Delete[tetatempprim, i ], i-1 ->  tetatnewtemp];
zcentertemp3 =ReplacePart[ Delete[zcenter1, i ], i-1 -> centernew];},
{NuQ1=Append[NuQtemp,0];
MuQ1=Append[MuQtemp,0];
l1temp4 = Append[l1temp, lnew];
l2temp4 =Append[l2temp,lnew/l2ratio];
tetatemp4 =Append[tetatemp,tetatnew];
tetatempprim4 = Append[tetatempprim, tetatnewtemp];
zcentertemp3 =Append[zcenter1, centernew];}]];

s5 = OpenWrite[OutputFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot+sign];
  Write[s5, Nu0];
  Write[s5, Mu0];
  Write[s5, s11];
  Write[s5, s12];
  Write[s5, s22];
  Write[s5, NuQ1];
  Write[s5, MuQ1];
  Write[s5, l1temp4];
  Write[s5, l2temp4];
  Write[s5, tetatemp4];
  Write[s5, tetatempprim4];
  Write[s5, zcentertemp3];
  Write[s5, zcentertemp3];
  Close[s5];]



(*delete inclusion i and put it at j place.*)
PlaceiToj[i_,j_] := {Print["--------------------Place Inclusion  ",i,"  to  ",j," ----------------------------------"];
Get["geometry.m"];
Needs["geometry`"];

tetatemp3=Delete[Insert[tetatemp,tetatemp[[i]],j+1],{{i}}];
tetatemp3prim=Delete[Insert[tetatempprim,tetatempprim[[i]],j+1],{{i}}];
l1temp3=Delete[Insert[l1temp,l1temp[[i]],j+1],{{i}}];
l2temp3=Delete[Insert[l2temp,l2temp[[i]],j+1],{{i}}];
NuQ3=Delete[Insert[NuQtemp,NuQtemp[[i]],j+1],{{i}}];
MuQ3=Delete[Insert[MuQtemp,MuQtemp[[i]],j+1],{{i}}];
zcentertemp3=Delete[Insert[zcenter1,zcenter1[[i]],j+1],{{i}}];
zcenternewtemp3=Delete[Insert[zcenter1,zcenter1[[i]],j+1],{{i}}];

  s5 = OpenWrite[OutputFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot];
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
  Close[s5];
Get["geometry.m"];
Needs["geometry`"];}



swap[a_,i_,j_]:=(ReplacePart[a,{i->a[[j]],j->a[[i]]}])


(*Replace1: When we have a crack penetration into a soft inclusion, the 
path of crack inside the inclusion is not a matter of interest for us. 
What we need at this point is if a crack emanates from the inclusion 
anymore? If yes, then we have to add more miro-crack to the system 
emanating from the inclusion. This means that I have to include soft 
inclusion as part of the crack path. This substitution can be done 
using the Replace1 function. *)

Replace1[i_] := {Print["---------------------Replace ",i," --------------------------------"];
Get["geometry.m"];
Needs["geometry`"];

Clear[tetatemp3,tetatemp3prim,l1temp3,l2temp3,NuQ3,MuQ3,zcentertemp3, zcenternewtemp3];
tetatemp3=Delete[ReplacePart[tetatemp,ntot->tetatemp[[i]]],{{i}}];
tetatemp3prim=Delete[ReplacePart[tetatempprim,ntot->tetatemp[[i]]],{{i}}];
l1temp3=Delete[ReplacePart[l1temp,ntot->l1temp[[i]]],{{i}}];
l2temp3=Delete[ReplacePart[l2temp,ntot->l2temp[[i]]],{{i}}];
NuQ3=Delete[ReplacePart[NuQtemp,ntot->NuQtemp[[i]]],{{i}}];
MuQ3=Delete[ReplacePart[MuQtemp,ntot->MuQtemp[[i]]],{{i}}];
zcentertemp3=Delete[ReplacePart[zcenter1,ntot->zcenter1[[i]]],{{i}}];
zcenternewtemp3=Delete[ReplacePart[zcenter1,ntot->zcenter1[[i]]],{{i}}];


  s5 = OpenWrite[OutputFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot -1];
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




Replaceij[i_,j_] := {Print["-------------------Replace ",i," and ", j," --------------------------------"];
Get["geometry.m"];
Needs["geometry`"];


tetatemp3=swap[tetatemp,i,j];
tetatemp3prim=swap[tetatempprim,i,j];
l1temp3=swap[l1temp,i,j];
l2temp3=swap[l2temp,i,j];
NuQ3=swap[NuQtemp,i,j];
MuQ3=swap[MuQtemp,i,j];
zcentertemp3=swap[zcenter1,i,j];
zcenternewtemp3=swap[zcenter1,i,j];

  s5 = OpenWrite[OutputFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot ];
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
  Close[s5];
Get["geometry.m"];
Needs["geometry`"];}


ReplaceInc[j_, ttest_] := 
 Block[{f}, f = For[i = 1, i < ttest + 1, i++,
    Get["geometry.m"];
    Needs["geometry`"];
    Replaceij[j, ntot - i]]];


DeleteInc[i_] := {Print["-------------------Delete Inclusion",i," -------------------------------"];
Get["geometry.m"];
Needs["geometry`"];

tetatemp3=Delete[tetatemp,{{i}}];
tetatemp3prim=Delete[tetatempprim,{{i}}];
l1temp3=Delete[l1temp,{{i}}];
l2temp3=Delete[l2temp,{{i}}];
NuQ3=Delete[NuQtemp,{{i}}];
MuQ3=Delete[MuQtemp,{{i}}];
zcentertemp3=Delete[zcenter1,{{i}}];
zcenternewtemp3=Delete[zcenter1,{{i}}];



  s5 = OpenWrite[OutputFile];
  Write[s5, n];
  Write[s5, num1];
  Write[s5, num2];
  Write[s5, ntot -1];
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
  Close[s5];
Get["geometry.m"];
Needs["geometry`"];}



(*This deletes output, copy outputtemp to output, and then deletes outputtemp.*)
OutputtempFile = FileNameJoin[{FilePath, "output", "outputtemp"}];
OutputTemp[]:={
DeleteFile[OutputFile];
ClearSystemCache[];
CopyFile[OutputtempFile,OutputFile];
DeleteFile[OutputtempFile];};


Copyoutput:=CopyFile[OutputFile,OutputtempFile];



(*Function r is defined first to copy the output file from first round 
simulation, and then append the coordinate of new inclusion without changing 
their size into the output1 file. Using this function crack propagation path 
can be determined. It is not given an output that the simulation is running 
on, but it will simply determined the new particle's coordinates that are 
added to the initial configuration. 
*)
(*inc is the number of inclusions left from the end of the crack that are not consider before merging.*)
routput1[t_,inc_,i_]:={(*Print["--------------------------routput1 i ----------------------------------"];*)
Get["geometry.m"];
Needs["geometry`"];

Output1File = FileNameJoin[{FilePath, "output", "output1"}];
If[t==1,CopyFile[OutputFile,Output1File],
Print["t is :  ",t];
s3=OpenRead[Output1File];
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

(*NuQ2=Join[Drop[NuQtemp1,-1-inc],Take[NuQtemp,-(inc+2)]];
MuQ2=Join[Drop[MuQtemp1,-1-inc],Take[MuQtemp,-(inc+2)]];
l1temp2=Join[Drop[l1temp1,-1-inc],Take[l1temp,-(inc+2)]];
l2temp2=Join[Drop[l2temp1,-1-inc],Take[l2temp,-(inc+2)]];
tetatemp2=Join[Drop[tetatemp1,-1-inc],Take[tetatemp,-(inc+2)]];
tetatemp2prim=Join[Drop[tetatemp1prim,-1-inc],Take[tetatempprim,-(inc+2)]];
zcentertemp2=Join[Drop[zcenter11,-1-inc],Take[zcenter1,-(inc+2)]];
zcenternewtemp2=Join[Drop[zcenter11,-1-inc],Take[zcenter1,-(inc+2)]];*)

NuQ2=Join[Drop[NuQtemp1,0-inc],Take[NuQtemp,-(inc+i)]];
MuQ2=Join[Drop[MuQtemp1,0-inc],Take[MuQtemp,-(inc+i)]];
l1temp2=Join[Drop[l1temp1,0-inc],Take[l1temp,-(inc+i)]];
l2temp2=Join[Drop[l2temp1,0-inc],Take[l2temp,-(inc+i)]];
tetatemp2=Join[Drop[tetatemp1,0-inc],Take[tetatemp,-(inc+i)]];
tetatemp2prim=Join[Drop[tetatemp1prim,-0-inc],Take[tetatempprim,-(inc+i)]];
zcentertemp2=Join[Drop[zcenter11,-0-inc],Take[zcenter1,-(inc+i)]];
zcenternewtemp2=Join[Drop[zcenter11,-0-inc],Take[zcenter1,-(inc+i)]];
DeleteFile[Output1File];
s4=OpenWrite[Output1File];

Write[s4,n];
Write[s4,num1];
Write[s4,num2];
Write[s4,ntemp+i];
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


(* This function reads data from the output file and delete all involved 
inclusions in the crack propagation path except the first two of them. 
It copies output to outputtemp, delete output, and then writes the new 
list to the output file. 
The aim of this part is to recalculate stress fields and SIFs for better 
accuracy using "naccuracy" variable to define the desire number of 
harmonic term.*)

Length2Crack[naccuracy_]:={Print["-------------------------Length2Crack ----------------------------------"];
s6=OpenRead[OutputFile];
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

CopyFile[OutputFile,OutputtempFile];
DeleteFile[OutputFile];
s7=OpenWrite[OutputFile];

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

(*SIF of pth inclusion caused by the qth one.*)
ExportSIF2[p_,q_]:= Block[{f,f1,kk,g1,g2,i,fg,k2}, 
Print["ExportSIF2[p_,q_]"];
PackageReload["run1`", KillShadowing -> True];
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];



f1 = Sum[l1temp[[i]], {i, ninitial, Times @@ Dimensions[l1temp]}];
g1=GQ[p];
Print["gc = ", g1]; 
kk = SIF`kcrack[p,q];
fg = (Chi0 + 1)/8/Mu0;(*1/GPa*)
k2 = kk*Conjugate[kk];
  g2 = Re[1000*fg*k2];Print["g[ntot] = ", g2];
PutAppend[List[Re[kk],Im[kk],g1,g2,f1],crackenergyOutputFile];
  f = If[g2 < g1, Print["Crack is arrested!"];, 
    Print["Crack is propagating!"]]];



End[] 
Protect@@Names["Changesize`*"]
EndPackage[] 

BeginPackage["Changesize`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","CombinMMBE`","eta0p`","EnergyTest`","SIF`"}] 
EndPackage[] 





















































