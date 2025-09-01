(* ::Package:: *)

BeginPackage["MergeFinal`"]
Unprotect@@Names["MergeFinal`*"];
ClearAll@@Names["MergeFinal`*"];




(*r3[]
The "r3" function writes all old and new inclusions of the system in the output file; update output file. 
*)


MergeFinal::usage="";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["MergeInclusion.m"]
Needs["MergeInclusion`"]




r3[]:={s2=OpenWrite[OutputFile];

Write[s2,n];
Write[s2,num1];
Write[s2,num2];
Write[s2,ntot1];

Write[s2,Nu0];
Write[s2,Mu0];

Write[s2,s11];
Write[s2,s12];
Write[s2,s22];

Write[s2,NuQ1];
Write[s2,MuQ1];

Write[s2,l1temp1];
Write[s2,l2temp1];
Write[s2,tetatemp1 ];
Write[s2,tetatemp1prim];
Write[s2,zcentertemp];
Write[s2,zcenternewtemp1];
Close[s2];
}


Needs["SubKernels`LocalKernels`"]
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]])
initsub[]

(*Merge inclusion p and q; to merge ntot-1 and ntot-2
MergeInclusion[ntot-1,ntot-2]*)

MergeFinal[p1_,q1_]:={

If[Abs[tetatemp[[q1]]- tetatemp[[p1]]] >= 10.0,Print["No Merge!"] ;Break[],
Print["Merge inclusions for " ,1 ," times"];
p=Min[p1,q1];
q=Max[p1,q1];
zc=zcenter[q]+l1[q] Exp[I*Thetaprim[q]];
za=zcenter[p]-l1[p] Exp[I*Thetaprim[p]];
zb=zcenter[q]-l1[q] Exp[I*Thetaprim[q]];
zbt=zcenter[p]+l1[p] Exp[I*Thetaprim[p]];
zbtt=zb-zbt;
rca=zc-za;
rc=Abs[rca]/2;
thetacarc=myArcTan[Re[rca],Im[rca]];
thetac=thetacarc*180/Pi;
thetac1 =If[thetac>90.0,thetac-180,If[thetac<-90,180+thetac,thetac]];
centernew=zc-rca/2.0;

Print["za=   " , za];
Print["zc=  " , zc];
Print["Delta= " , zbtt];

tetatemp1 =Insert[Delete[tetatemp,{{q},{p}}],thetac1,+Min[p,q]];
tetatemp1prim =Insert[Delete[tetatempprim,{{q},{p}}],thetac,+Min[p,q]];
l1temp1 = Insert[Delete[l1temp,{{q},{p}}],rc,+Min[p,q]];
l2temp1 =Insert[Delete[l2temp,{{q},{p}}],rc/l2ratio,+Min[p,q]];
NuQ1 =Insert[Delete[NuQtemp,{{q},{p}}],0,+Min[p,q]];
MuQ1 =Insert[Delete[MuQtemp,{{q},{p}}],0,+Min[p,q]];
zcentertemp = Insert[Delete[zcenter1,{{q},{ntot-1}}],centernew,+Min[p,q]];
zcenternewtemp1 = Insert[Delete[zcenternewlist1,{{q},{p}}],centernew,+Min[p,q]];
ntot1=ntot-1;
r3[];]}


End[] 
Protect@@Names["MergeFinal`*"]
EndPackage[] 

BeginPackage["MergeFinal`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","InterfaceFunc`","EllipseInterfaceCheck`","EnergyTest`","LineInterfaceCheck`","Addone`","AddInclusion`","MergeInclusion`"}] 
EndPackage[] 











































