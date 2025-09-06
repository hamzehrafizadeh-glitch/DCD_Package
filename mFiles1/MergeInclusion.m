(* ::Package:: *)

BeginPackage["MergeInclusion`"]
Unprotect@@Names["MergeInclusion`*"];
ClearAll@@Names["MergeInclusion`*"];




(*r3[]
The "r3" function writes all old and new inclusions of the system in the output file; update output file. 

Main Functions
2-MergeInclusion[]
The "MergeInclusion" function merges two inclusions "ntot-1", and "ntot-2" in the system if their relative inclination angle is less than a predetermined amount. The merging process is stopped by either not satisfying the relative inclination angle condition or by elongating the crack length to more than "?" times of the total length of the original crack. 
*)


MergeInclusion::usage="function merges two inclusions ntot-1, and ntot-2 in the system if their relative inclination angle is less than a predetermined amount. The merging process is stopped by either not satisfying the relative inclination angle condition or by elongating the crack length to more than ? times of the total length of the original crack. ";



Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["AddInclusion.m"]
Needs["AddInclusion`"]




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


(*Needs["SubKernels`LocalKernels`"]
initsub[]:=(subkernel=LaunchKernels[LocalMachine[4]])
initsub[]*)

(*Merge inclusion p and q; to merge ntot-1 and ntot-2
MergeInclusion[ntot-1,ntot-2]*)

MergeInclusion[p1_,q1_]:={
j=1;i=0;While[j==1,i++;
p=Min[p1,q1];
q=Max[p1,q1];
If[Abs[tetatempprim[[q1]]- tetatempprim[[p1]]] >=0.5,Print["Delta Theta is: ",Abs[tetatempprim[[q1]]- tetatempprim[[p1]]]," No Merge!"] ;Break[],
Print["Merge inclusions for " ,i ," times"];
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
If[Abs[zbtt]>1.1delta, Print[Style["Delta is Larger than 0.2, it cause an uncertainty in the results!", FontColor -> Purple]];Break[]];

tetatemp1 =Insert[Delete[tetatemp,{{q},{p}}],thetac1,+Min[p,q]];
tetatemp1prim =Insert[Delete[tetatempprim,{{q},{p}}],thetac,+Min[p,q]];
l1temp1 = Insert[Delete[l1temp,{{q},{p}}],rc,+Min[p,q]];
l2temp1 =Insert[Delete[l2temp,{{q},{p}}],rc/100,+Min[p,q]];
NuQ1 =Insert[Delete[NuQtemp,{{q},{p}}],0,+Min[p,q]];
MuQ1 =Insert[Delete[MuQtemp,{{q},{p}}],0,+Min[p,q]];
zcentertemp = Insert[Delete[zcenter1,{{q},{p}}],centernew,+Min[p,q]];
zcenternewtemp1 = Insert[Delete[zcenternewlist1,{{q},{p}}],centernew,+Min[p,q]];
ntot1=ntot-1;
r3[];

AddInclusion[1];

lref=Sum[l1temp[[i]],{i,ntot}];
Print["lrefAdd:   ",lref];
If[(lref/ltemp)>=3.0,Print["Crack propagates ", 2," times as much as its previous length!"];Break[]]
]]}


End[] 
Protect@@Names["MergeInclusion`*"]
EndPackage[] 

BeginPackage["MergeInclusion`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","InterfaceFunc`","EllipseInterfaceCheck`","EnergyTest`","LineInterfaceCheck`","Addone`","AddInclusion`"}] 
EndPackage[] 











































