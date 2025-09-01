(* ::Package:: *)

BeginPackage["run2`"]
Unprotect@@Names["run2`*"];
ClearAll@@Names["run2`*"];


(*Note
There is a problem in our calculation especially when the main crack crosses the 
other inclusions. The system encounter problem as the main crack pass the branch cut. 
This problem might be revealed by replacing or omitting the inclusion with the 
involved interface. 

I have to check this problem later. This might be easily solved by omitting the 
crossed inclusion; since the cross inclusion should not have significant effect 
except shadowing in this case. I can find the crossed inclusion by checking if 
the crack has any intersect with the branch cuts. *)


(*main-crack-One Directional propagation RHS *)


(*This file produces the propagation path direction for the main crack. 
It encompasses five functions in which two of them are the main functions. 

1- eta0[p]
The "eta0" function finds a propagation angle from the original inclusion p, using the maximum tangential stress criteria.

2-rtest[p]
The "rtest" function determines the direction of the new inclusion in the global coordinate system using the fact that the crack propagation is in the perpendicular direction to the surface of the original crack/inclusion. The "lnew" is the distance between the new inclusion center and the surface of the original inclusion/crack. 

3-r3[]
The "r3" function writes all old and new inclusions of the system in the output file; update output file. 

Main Functions

1-AddInclusion[t]
The "AddInclusion" function adds "t" number of inclusions with a specific length "lnew1" to the system. Adding process is started from the last inclusion "ntot-th" in the input file. 

2-MergeInclusion[]
The "MergeInclusion" function merges two inclusions "ntot-1", and "ntot-2" in the system if their relative inclination angle is less than a predetermined amount. The merging process is stopped by either not satisfying the relative inclination angle condition or by elongating the crack length to more than "?" times of the total length of the original crack. 
*)


(*eta0::usage="It determines an angle, in which the tangential stress on the surface of the inclusion q is maximum. "; *)
lnew1::usage="lnew1/2 is a vertical distance from the surface of the inclusion at eta0 angle."; 
delta::usage="It defines the distance between new added inclusions.";
AddInclusion::usage="This function adds t number of inclusions with a specific length lnew1 to the system. Adding process is started from the last inclusion ntot-th in the input file.";
MergeInclusion::usage="function merges two inclusions ntot-1, and ntot-2 in the system if their relative inclination angle is less than a predetermined amount. The merging process is stopped by either not satisfying the relative inclination angle condition or by elongating the crack length to more than ? times of the total length of the original crack. ";
NewCenter::usage="Center of new inclusion to add to the system as crack propagates.";
AddOne::usage="This function adds 1 inclusion with a specific length lnew1 to the system. Adding process is started from the last inclusion ntot-th in the input file.";
MergeFinal::usage="";
g::usage="Maximum Energy release rate of a propagating crack.:";
gc::usage="Maximum Fracture energy at the crack tipe..";
PropagationTest::usage="Check if the crack energetically can propagate?:";
myArcTan::usage="myArcTan[x_Real], it returns the angle in a range of [0,2Pi]. ";
Replace2::usage="";
etadeltaGB::usage="";
l2ratio::usage="It is the 1/e ratio, where e is the aspect ratio.";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];
*)


Get["PackageManipulations.m"]
Needs["PackageManipulations`"]
Get["run1.m"]
Needs["run1`"];
Get["SIF.m"]
Needs["SIF`"];



myArcTan[x_,y_] := Piecewise[{{Pi + ArcTan[y/x], x < 0}}, ArcTan[y/x]];


eta0[p_]:=Block[{r, r1, r2,r3},Unprotect[eta];Clear[Xi];
 r1 = SigmaEtaEta[p];
 Xi = 1.*Zeta0[p] + I*eta;
 r2 = SEtaEta[p]; 
 r3 = r1 + r2/.{Conjugate[eta] -> eta};
 r=If[BE==1,Re[r3],Re[r1]];
f1= If[0 <= Abs[Thetaprim[p]] <= Pi/2, 
    Last[NMaximize[{r, -Pi/2 < eta < Pi/2}, eta, AccuracyGoal -> 5]], 
     Last[NMaximize[{r, Pi/2 < eta < 3*Pi/2}, eta, AccuracyGoal -> 5]]];
f=eta /.f1]

(*By selecting r2 as a propagation angle,the crack propagation is confined to merely 
the RHS propagation.This means that the angle is between [-Pi/2 ,Pi/2]*)

(*In order to insure that the energy release associated with the crack extension is 
maximum,the sign of the propagation angle $\theta$ should be opposite to the sign of
$K_{II}$ \cite{wei1982nonlinear}*)

eta01[p_]:=Block[{f,r1,r2,k11,k22},If[ntot==1,q=p,q=p-1];
k11=Re[kcrack[p,q]];
k22=Im[kcrack[p,q]];
r1=ArcCos[(3*k22^2-Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
r2=ArcCos[(3*k22^2+Sqrt[k11^4+8.0*(k11*k22)^2])/(k11^2+9*k22^2)];
f=If[k22>0,-r2,r2]]
(*For this one I have to write a new NewCenter[p_] one considering our new angle*)





(*It obtaines the position of the new coordinate system origin where the new inclusion is
originated from the pth inclusion.*)

rtest[p_]:=Block[{f,d},Clear[x,y];d=Re[dd[p]*Exp[-I Theta[p]]];
f={x,y}/.Solve[
{d*(Cosh[Zeta0[p]]^2-Sinh[Zeta0[p]]^2) Cos[feta] Sin[feta]
-x*Cosh[Zeta0[p]]*Sin[feta]+y*Sinh[Zeta0[p]]*Cos[feta]==0,
((x-d* Cosh[Zeta0[p]] Cos[feta])^2+
(y-d* Sinh[Zeta0[p]] Sin[feta])^2)==lnew1^2,(x/(d*Cosh[Zeta0[p]]))^2
+(y/(d* Sinh[Zeta0[p]]))^2-1>0},{x,y}]]

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
(*initsub[]*)




(*This is a function to add t new inclusions after ntot inclusion.*)
ltemp=6.0;
Print["ltemp ", ltemp];
lnew1=0.6;
delta = 0.2;
lnew=lnew1-delta;
Print["lnew ",lnew];
l2ratio=100;

NewCenter[p_]:=Block[{f1,r,r1,f,xx1,yy1,x0,y0},r=rtest[p];
x0=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[feta];
y0=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[feta];

If[Dimensions[r]=={1,2},xx1=r[[1]][[1]];yy1=r[[1]][[2]],Table[s=If[Sign[r[[i]][[1]]]==Sign[x0]&&Sign[r[[i]][[2]]]==Sign[y0],1,0];
If[s==1,xx1=r[[i]][[1]];yy1=r[[i]][[2]]],{i,2}]];
f=xx1+I*yy1]
(*NewCenter[ntot] is equal to (xx1+I*yy1)*)


(*Addone is a function that Uses eta0,and NewCenter to add one inclusion to the 
system. The concept is based on the maximum shear stress and using SIFs to find a
propagation direction.*)
AddOne[p_]:=Block[{r,x0,y0,xx1,yy1},feta=eta0[p];
x0=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[feta];
y0=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[feta];

r=NewCenter[p];
xx1=Re[r];
yy1=Im[r];

rt1=Thetaprim[p]*180/Pi+Re[N[ArcTan[(yy1-y0)/(xx1-x0)]/\[Degree]]];
rt=If[0.00001<=rt1<=270 || -90<=rt1<=0 ,rt1,rt1-360];
Print["rt is an inclination angle: ", rt];
Thetanew =If[rt>90.0001,rt-180,If[rt<-90.01,180+rt,rt]];
Print["Thetanew is: ",Thetanew];
Thetaarc=Thetanew*Pi/180.0;
zcenternew=r*Exp[I Theta[p]]+zcenternewlist1[[p]];

    tetatemp1 = Append[tetatemp, Thetanew];
    tetatemp1prim = Append[tetatempprim, rt];
    l1temp1 = Append[l1temp, lnew];
    l2temp1 = Append[l2temp, lnew/l2ratio];
    NuQ1 = Append[NuQtemp, 0.0];
    MuQ1 = Append[MuQtemp, 0.0];
    zcentertemp = Append[zcenter1, zcenternew];
    zcenternewtemp1 = Append[zcenternewlist1, zcenternew]; 
    ntot1 = ntot + 1;r3[]]






(*AddInclusion[t_]:={For[q = 1, q < t + 1, q++, Print["AddInclusion number: ",q];
	 
CloseKernels[subkernel]; initsub[];  
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];

AddOne[ntot];
If[q==t,Get["geometry.m"];Needs["geometry`"]];]}*)


(*Replace1: When we have a crack penetration into a soft inclusion, the 
path of crack inside the inclusion is not a matter of interest for us. 
What we need at this point is if a crack emanates from the inclusion 
anymore? If yes, then we have to add more miro-crack to the system 
emanating from the inclusion. This means that I have to include soft 
inclusion as part of the crack path. This substitution can be done 
using the Replace1 function. *)

Replace1[i_] := {Print["------------------------Replace ",i," ----------------------------------"];
Get["geometry.m"];
Needs["geometry`"];


tetatemp3=Delete[ReplacePart[tetatemp,ntot->tetatemp[[i]]],{{i}}];
tetatemp3prim=Delete[ReplacePart[tetatempprim,ntot->tetatemp[[i]]],{{i}}];
l1temp3=Delete[ReplacePart[l1temp,ntot->l1temp[[i]]],{{i}}];
l2temp3=Delete[ReplacePart[l2temp,ntot->l2temp[[i]]],{{i}}];
NuQ3=Delete[ReplacePart[NuQtemp,ntot->NuQtemp[[i]]],{{i}}];
MuQ3=Delete[ReplacePart[MuQtemp,ntot->MuQtemp[[i]]],{{i}}];
zcentertemp3=Delete[ReplacePart[zcenter1,ntot->zcenter1[[i]]],{{i}}];
zcenternewtemp3=Delete[ReplacePart[zcenter1,ntot->zcenter1[[i]]],{{i}}];


  s5 = OpenWrite[
    OutputFile];
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





Replace2 := {Print["------------------------Replace ",2," ----------------------------------"];
Get["geometry.m"];
Needs["geometry`"];

tetatemp3 = ReplacePart[tetatemp,ntot->Thetanewp];
    tetatemp3prim = ReplacePart[tetatempprim,ntot->rtp];
   l1temp3=l1temp;
l2temp3=l2temp;
NuQ3=NuQtemp;
MuQ3=MuQtemp;
    zcentertemp3 = ReplacePart[zcenter1,ntot->centerprime];
    zcenternewtemp3 =zcentertemp3;

  s5 = OpenWrite[
    OutputFile];
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
  Close[s5];}





Clear[GrainB]
acorner=1;bcorner=10;xc=-2;yc=0;
GrainB[x_,y_]:=Block[{f},f=Abs[(x-xc)/acorner+(y-yc)/bcorner]+Abs[(x-xc)/acorner-(y-yc)/bcorner]]

GrainB[x,y]==2

GBIntersect[p_]:=Block[{f,f1,f2,f3,f4,f5,x1,y1,pluspoint,midpoint,minuspoint,m,cline,it,xx,yy},pluspoint=zcenter1[[p]]-l1[p]*Exp[I*Thetaprim[p]];
midpoint=zcenter1[[p]];
minuspoint=zcenter1[[p]]+l1[p]*Exp[I*Thetaprim[p]];
x1=Re[midpoint];
y1=Im[midpoint];
m=Tan[Thetaprim[p]];
cline=y1-m*x1;
f={yy==m*xx+cline,GrainB[xx,yy]==2};
f1=Solve[f,{xx,yy}];
f2=If[f1!={},{xx,yy}/.f1,{}];
If[f1!={},it=Times@@Dimensions[First[f2]];
f3=Sum[f2[[i]],{i,it}]/it;
f4=Table[If[GrainB[Re[i],Im[i]]<=2,1,0],{i,{pluspoint,midpoint,minuspoint}}];
f5=f4[[1]]+f4[[2]]+f4[[3]],f5=0];
If[f5==0,gb=ntot+2,Print["Grain Boundary!"];gb=1;f3]]

etadeltaGB[p_]:=Block[{x,y,fxy,fpxy,r,ftheta,rtp1,rtp2},
fxy=Sqrt[((x-xc)/acorner+(y-yc)/bcorner)^2]+Sqrt[((x-xc)/acorner-(y-yc)/bcorner)^2];
fpxy=\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]fxy\) +I*\!\(
\*SubscriptBox[\(\[PartialD]\), \(y\)]fxy\) ;
r=GBIntersect[p];

x=First[r];
y=Last[r];
ftheta=myArcTan[Im[fpxy],Re[fpxy]]*180/Pi;

rtp1=-90+ftheta;
rtp2=+90+ftheta;
rtp =If[Abs[rtp1- tetatempprim[[ntot-1]]]<=90,rtp1,rtp2];

Print["rtp is an inclination angle of Replace2: ", rtp];
Thetanewp =If[rtp>90.0,rtp-180,If[rtp<-90,180+rtp,rtp]];
Print["Thetanewp of Replace 2 is: ",Thetanewp];
Thetaarcp=Thetanewp*Pi/180.0 ;
centerprime=x+I*y+lnew1*Exp[I*(Pi*rtp/180)]]



etadelta[p_]:={PnewZ=zcenter1[[ntot-1]]+(l1[ntot-1])*Exp[I*Thetaprim[ntot-1]];
Pnew=(PnewZ-zcenter[p])*Exp[-I*Theta[p]];
d=Re[dd[p]*Exp[-I Theta[p]]];

xetadelta=Re[Pnew]* Cos[Theta[p]]-Im[Pnew]* Sin[Theta[p]];
yetadelta= Im[Pnew]* Cos[Theta[p]]+ Re[Pnew]* Sin[Theta[p]];
f={d*(Cosh[Zeta0[p]]^2-Sinh[Zeta0[p]]^2) Cos[etap] Sin[etap]
-Re[Pnew]*Cosh[Zeta0[p]]*Sin[etap]+Im[Pnew]*Sinh[Zeta0[p]]*Cos[etap]==0,
((Re[Pnew]-d* Cosh[Zeta0[p]] Cos[etap])^2+
(Im[Pnew]-d* Sinh[Zeta0[p]] Sin[etap])^2)==lnewp^2,(Re[Pnew]/(d*Cosh[Zeta0[p]]))^2
+(Im[Pnew]/(d* Sinh[Zeta0[p]]))^2-1>0};

f1=Solve[f,{etap,lnewp}];
f2={etap,lnewp}/.f1;
f3=Select[f2,Element[#,Reals]&];

idom=Times@@First[Dimensions[f3]];
jdom=Times@@Last[Dimensions[f3]];

expr=Table[If[f3[[i,j]]>0,f3[[i]]],{i,idom},{j,2,jdom}];
expr1=expr/.{Null}->Sequence[];
expr2=Replace[expr1,{x_List}:>x,{0,-3}];
yt=Max[Table[expr2[[i,2]],{i,Times@@First[Dimensions[expr2]]}]];
exprt=Table[If[expr2[[i,2]]!=yt,expr2[[i]]],{i,Times@@First[Dimensions[expr2]]}];
exprt1=exprt/.Null->Sequence[];
If[Dimensions[exprt1]=={2,2}&&exprt1[[1,2]]==exprt1[[2,2]],etanew=Abs[exprt1[[1,1]]];deltanew=exprt1[[1,2]],etanew=exprt1[[1,1]];
deltanew=exprt1[[1,2]]];

x0p=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[etanew];
y0p=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[etanew];


xx1p=Re[Pnew];
yy1p=Im[Pnew];

rtp1=-90+Re[N[ArcTan[(yy1p-y0p)/(xx1p-x0p)]/\[Degree]]]+Theta[p]*180/Pi;
rtp2=+90+Re[N[ArcTan[(yy1p-y0p)/(xx1p-x0p)]/\[Degree]]]+Theta[p]*180/Pi;
rtp =If[Abs[rtp1- tetatempprim[[ntot-1]]]<=90,rtp1,rtp2];(* If[Abs[rtp1]<=90,rtp1,rtp2];*)

Print["rtp is an inclination angle of Replace2: ", rtp];
Thetanewp =If[rtp>90.0,rtp-180,If[rtp<-90,180+rtp,rtp]];
Print["Thetanewp of Replace 2 is: ",Thetanewp];
Thetaarcp=Thetanewp*Pi/180.0 ;


dist=0.1;
m=Tan[Re[N[ArcTan[(yy1p-y0p)/(xx1p-x0p)]/\[Degree]]]+Theta[p]*180/Pi];
c=Im[PnewZ]-m*Re[PnewZ];
delttemp=deltanew-dist;
deltap=If[Abs[delttemp]<=0.001,0.001 ,delttemp];

f1=Solve[{(Im[PnewZ]-yy)^2 +(Re[PnewZ]-xx)^2 ==deltap^2,
yy==m*xx+c},{xx,yy}];
f2={xx,yy}/.f1;

gtemp=d* Cosh[Zeta0[p]] Cos[etanew]+I*(d* Sinh[Zeta0[p]] Sin[etanew])+zcenter[1];

q1=((f2[[1,1]]-Re[gtemp])^2+
(f2[[1,2]]-Im[gtemp])^2);

q2=((f2[[2,1]]-Re[gtemp])^2+
(f2[[2,2]]-Im[gtemp])^2);

cent=If[q1<q2 ,f2[[1,1]]+I*f2[[1,2]],f2[[2,1]]+I*f2[[2,2]]];

centerprime=cent+lnew1*Exp[I*(Pi*rtp/180)]*Exp[0*I Theta[p]];}



readstart[]:=Block[{s8,rtest},s8 = OpenRead[
   FileNameJoin[{FilePath, "output", "start"}]];
rtest =Flatten[ Table[ReadList[s8, Expression, 1], {j, 1+1}]]]

ninitial=readstart[][[1]];
lltest=readstart[][[2]];
(*I have to find a way for avoiding or adding ttest in the code. ttest is used in the run4 file.*)


(*Note: inc<=ttest-2*)(*decide how many inclusions left at the end of the crack path to 
predict path more accurate.*)



(*I have to check intersection with all initial inclusions except the one which is the initial crack*)
(*r=Table[If[Testintersect[ntot,i]\[NotEqual]0,Print[ntot, " and ",i, " inclusions intersect!"];2,1],{i,ntot-1}];
*)

rreplace:=Block[{r,f},r=For[i=1,i<(ninitial),i++,
     
If[Testintersect[ntot,i]!=0,Print[ntot, " and ",i, " inclusions intersect!"];Return[i];Break[],Return[ntot+2]]]; 
f=If[ninitial==1,Return[ntot+2];Print["yes"],r]]



(*gGB=0.9*gc;
AddInclusion[t_]:={For[q=1,q<t+1,q++,Print["AddInclusion number: ",q];
j=1;While[j==1,
PackageReload["run1`",KillShadowing->True];
CloseKernels[subkernel];initsub[];
PackageReload["SIF`",KillShadowing->True];
CloseKernels[subkernel];initsub[];
AddOne[ntot];
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
GBIntersect[ntot];
rtemp=rreplace;(*rreplace;*)
If[rtemp!=(ntot+2),If[(MuQ[rtemp]/Mu0)>1.,etadelta[rtemp];Replace2;j=2,Replace1[rtemp]],If[gb!=(ntot+2) && gGB<gc,etadeltaGB[ntot];Replace2;j=2,j=2]]];]}
*)
gGB=0.9*gc;
AddInclusion[t_]:={For[q=1,q<t+1,q++,Print["AddInclusion number: ",q];
j=1;While[j==1,
PackageReload["run1`",KillShadowing->True];
CloseKernels[subkernel];initsub[];
PackageReload["SIF`",KillShadowing->True];
CloseKernels[subkernel];initsub[];
AddOne[ntot];
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
GBIntersect[ntot];
rtemp=rreplace;(*rreplace;*)

r1=If[rtemp==(ntot+2),0,1];
r2=If[gb==(ntot+2),0,1];


If[r1==0 && r2==0,j=2,
If[r1==1,If[r2==1,If[gGB < gc,etadeltaGB[ntot];Replace2;j=2,If[(MuQ[rtemp]/Mu0)>1.,etadelta[rtemp];Replace2;j=2,Replace1[rtemp]]],If[(MuQ[rtemp]/Mu0)>1.,etadelta[rtemp];Replace2;j=2,Replace1[rtemp]]],If[gGB < gc,etadeltaGB[ntot];Replace2;j=2,j=2]]
]

];]}


(*here I imagine that in the system under calculation,a crack path consists of two micro-cracks.*)

(*At this point and before determining whether the crack can propagate,deflect,or arrest,
it should be determined where is the crack tip and what is the corresponding fracture energy.*)

Iinterface:=Block[{r,f},r=For[i=1,i<ninitial,i++,     
If[Testinterface[ntot,i]!=0,(*Print[ntot, " is in the ",i, "th inclusions interface!"];*)Return[i];Break[],Return[ntot+2]]];
f=If[ninitial==1,Return[ntot+2];Print["yes"],r]]

g[p_, ttest_] := 
  Block[{fg, kk, k2, f}, fg = (Chi0 + 1)/8/Mu0;(*1/GPa*)
   kk = If[ntot == ninitial, SIF`kcrack[p, p], 
     SIF`kcrack[p - ttest, p - ttest + 1]];
   k2 = kk*Conjugate[kk];
   f = Re[1000*fg*k2]];
(*At this point and before determining whether the crack can propagate,deflect,or arrest,it should be determined where is the crack tip and what is the corresponding fracture energy.*)

gc:=Block[{itemp},
Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
itemp = Iinterface;
(*Here I assumed that at the first step, crack is not close to the interface!*)
If[itemp!=(ntot+2),GintQ[itemp],G0]]

(*l1temp1[[i]] is defined in cracklength function in run4*)
PropagationTest[p_, ttest_] := 
 Block[{f,f1,kk,g1,g2}, f1 = Sum[l1temp[[i]], {i, ninitial, Times @@ Dimensions[l1temp]}];
g1=gc;g2=g[p, ttest];
Print["gc = ", g1]; Print["g[ntot] = ", g2];
kk = If[ntot == ninitial, SIF`kcrack[p, p], 
     SIF`kcrack[p - ttest, p - ttest + 1]];
PutAppend[List[Re[kk],Im[kk],g1,g2,f1],crackenergyOutputFile];
  f = If[g2 < g1, Print["Crack is arrested!"];, 
    Print["Crack is propagating!"];]];


AddInclusiontest[t_]:={For[q = 1, q < t + 1, q++, Print["AddInclusion number: ",q];
	 
j=1;While[j==1,
CloseKernels[subkernel]; initsub[];  
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];

rtemp=If[ntot==1,ntot+2,rreplace];
Print["rtemp = ", rtemp]
Print["gc = ", gc];Print["g[ntot] = ", g[ntot]];
If[g[ntot]<gc,Print["Crack is arrested!"];Abort[],Print["Crack is propagating!"];];
AddOne[ntot];

Get["expansionCoefficients.m"];
Needs["ExpansionCoefficient`"];
Clear[rtemp];
rtemp = If[ntot==1,ntot+2,rreplace];(*rreplace;*)
If[rtemp!=(ntot+2),If[(MuQ[rtemp]/Mu0)>1.,etadelta[rtemp];Replace2;j=2,Replace1[rtemp]],j=2]]];
}


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
Protect@@Names["run2`*"]
EndPackage[] 

BeginPackage["run2`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","PackageManipulations`","run1`"}] 

EndPackage[] 



























































































