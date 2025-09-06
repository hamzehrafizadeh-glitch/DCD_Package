(* ::Package:: *)

BeginPackage["EllipseInterfaceCheck`"]
Unprotect@@Names["EllipseInterfaceCheck`*"];
ClearAll@@Names["EllipseInterfaceCheck`*"];


(*Test Interface is not a symmetric function.the qth inclusion has 
an interface function. This function finds if two ellipses have an intersect.*)






(*Testintersect::usage=" Testinterface[p_,q_]: It is a function to check whether the inclusion p, overlap or intersect with the inclusion q and its interface?";*)
IntersectCheck::usage="It is a function to check whether any pair of inclusions in the system intersect? ";
Testinterface::usage="IntersectCheck: Testinterface[p_,q_]: It is a function to check whether the inclusion p, overlap or intersect with the interface of the inclusion q?";
EllipseIntersectPoint::usage="EllipseIntersectPoint[i_,p_,tt_]: This function calculates a line and an ellipse intersection point.";
tt::usage="tt[i_]: This value identifies the interface domain.";
Testinterface1::usage="Testinterface1[j_]:";
IntersectCheckQ::usage="IntersectCheckQ[q_]";
TestPenetration::usage="TestPenetration[j]";
PassInclusion::usage="PassInclusion[pp1_]";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Needs["eta0p`"];
Needs["Changesize`"];
Needs["InterfaceFunc`"];

(*-----------------------------INTERSECT CHECK---------------------------*)
(*Solve[x*dd[1]-dist\[Equal]lnew1*Cos[75Degree],{x}]*)


Clear[x1,y1,xxx,yyy]
x1[i_]:=(Re[zte]-Re[zcenter[i]])*Cos[Theta[i]]+(Im[zte]-Im[zcenter[i]])*Sin[Theta[i]];
y1[i_]:=(Im[zte]-Im[zcenter[i]])*Cos[Theta[i]]-(Re[zte]-Re[zcenter[i]])*Sin[Theta[i]];

xxx[i_,q_]:=Block[{f},f=(Re[i]-Re[zcenter[q]])*Cos[Theta[q]]+(Im[i]-Im[zcenter[q]])*Sin[Theta[q]]]
yyy[i_,q_]:=Block[{f},f=(Im[i]-Im[zcenter[q]])*Cos[Theta[q]]-(Re[i]-Re[zcenter[q]])*Sin[Theta[q]]]


(*(tt*dd)^2=l1^2+l2^2=(1+cof^2)l1prim^2 dd^2
=(1+cof^2)l1^2 (1+cof^2)=dd^2/l1^2 (tt*dd)^2
=(1+cof^2)l1prim^2=(dd^2/l1^2)l1prim^2 
Thus l1prim=tt*l1*)(*Test Intersect is not a symmetric function.the qth inclusion has 
an interface function.*)

(*Testintersect[p_,q_]:=Block[{f,f1,f2,f3,f4,f5,tt1,pluspoint,minuspoint},If[(Abs[zcenter[q]-zcenter[p]]>(l1[p]+tt*l1[q]))||p\[Equal]q,0,f=Table[tt1=If[i\[Equal]q,tt,1];((x1[i])/Cosh[Zeta0[i]])^2+((y1[i])/Sinh[Zeta0[i]])^2\[Equal]Re[(tt1*dd[i])*Conjugate[(tt1*dd[i])]],{i,{p,q}}];
f1=Solve[f,{zte}];
f2=If[f1\[NotEqual]{},{zte}/.f1,{}];
f3=First[Dimensions[f2]];
pluspoint=zcenter1[[p]]-l1[p]*Exp[I*Thetaprim[p]];
minuspoint=zcenter1[[p]]+l1[p]*Exp[I*Thetaprim[p]];
f4=Table[If[((xxx[i,q])/Cosh[Zeta0[q]])^2+((yyy[i,q])/Sinh[Zeta0[q]])^2-Re[(tt*dd[q])*Conjugate[(tt*dd[q])]]<0,1,0],{i,{pluspoint,minuspoint}}];
f5=f4[[1]]+f4[[2]]+f3]]*)

(*How to define tt: 
l1=Cosh[Zeta0]*Dq;
l2=Sinh[Zeta0]*Dq;
tt=1+?
If I want the distance to be equal to delta, then
delta\[Equal]?*l1
tt=1+delta/l1*)

Clear[tt]
tt[i_]:=1+1.*delta/l1[i];(*delta/l1[i]*)


Testinterface[p_,q_]:=Block[{zte,i,f,f1,f2,f3,f4,f04,f5,f6,tt1,ti,pluspoint,minuspoint,ltest0,test},ltest0 = 0.5*(2.0*lnew + delta);ti=tt[q];
If[(Abs[zcenter[q]-zcenter[p]]>(l1[p]+ti*l1[q]))||p==q,{0,0,0,0},f=Table[tt1=If[i==q,ti,1];((x1[i])/Cosh[Zeta0[i]])^2+((y1[i])/Sinh[Zeta0[i]])^2==Re[(tt1*dd[i])*Conjugate[(tt1*dd[i])]],{i,{p,q}}];
f1=Solve[f,{zte}];
f2=If[f1!={},{zte}/.f1,{}];
f3=First[Dimensions[f2]];
Print["Intersection points are  ", f2];
(*Print["f3  ",f3];*)
minuspoint=StartPoint[p];
pluspoint=EndPoint[p];
Print["pluspoint  ", pluspoint];
f4=Table[If[((xxx[i,q])/Cosh[Zeta0[q]])^2+((yyy[i,q])/Sinh[Zeta0[q]])^2-Re[(ti*dd[q])*Conjugate[(ti*dd[q])]]<=10^(-5),1,0],{i,{pluspoint,minuspoint}}];(*This shows tha the start and end point of the inclusion p is in the inclusion q.*)
f04=Table[test=((xxx[i,q])/Cosh[Zeta0[q]])^2+((yyy[i,q])/Sinh[Zeta0[q]])^2-Re[Dq[q]^2];If[Sign[test]==-1,1,0],{i,{pluspoint(*,minuspoint*)}}];

Print["HERE I CHANGED THE FUNCTION USING TEST!"];
f5=f4[[1]]+f4[[2]]+f3;
f6=f04[[1]](*+f4[[1]]*);
(*Print["f4  ", f4];*)
Print["{intersect,pluspoint,minuspoint,pluspoint inside the inclusion}   ", Flatten[{f5,f4,f6}]];
Flatten[{f5,f4,f6}]]
]


(*This function return the number of the first inclusion that intersect with another one. 
If no intersection is seen, then it is null! *)
(*IntersectCheck:=For[p=1,p\[LessEqual]ninitial,p++,For[q=p+1,q<=ntot,q++,Print["IntersectCheck"];
If[Testinterface[q,p]!=0,Print["Oops inclusions ",p," and ",q," intersect!"];Return[p],0]]]
*)
IntersectCheck:=Block[{f,f1,f2,f3,p},f1=For[p=1,p<=ninitial,p++,For[q=p+1,q<=ntot,q++,f2=Testinterface[q,p];f3=f2[[1]]*f2[[2]]*f2[[4]];
Print["f3   ",f3];
If[(f2[[1]]!=0 && f2[[2]]!=0),If[f3==0,Print["Oops inclusions ",p," and ",q," intersect!"];i={p};Break[],Print["Crack is permeating into the inclusion  ", p];i={p,3}];Break[],i={0}]]];
f=i]


IntersectCheckQ[q_]:=Block[{f,f1,f2,f3,i,p},If[ninitial==0,Print["ninitial =0 and No intersect!"];i={0},
f1=For[p=1,p<=ninitial,p++,f2=Testinterface[q,p];f3=f2[[1]]*f2[[2]]*f2[[4]];
Print["f2   ",f2];
If[(f2[[1]]!=0 && f2[[2]]!=0),If[f2[[4]]==0,Print["Oops inclusions ",p," and ",q," intersect!"];i={p};Break[],Print["Crack is permeating into the inclusion  ", p];i={p,3};Break[]],i={0}]]];
f=i]


EllipseIntersectPoint[i_,p_,tt_]:=Block[{x,y,r1,r2,Dis1,Dis2,x1,y1,t,t1,x0,y0,m,cline,f,f1,point},

x1=(x-Re[zcenter[i]])*Cos[Theta[i]]+(y-Im[zcenter[i]])*Sin[Theta[i]];
y1=(y-Im[zcenter[i]])*Cos[Theta[i]]-(x-Re[zcenter[i]])*Sin[Theta[i]];
t=EndPoint[p];
t1=StartPoint[p];
x0=Re[t];y0=Im[t];
m=Tan[Thetaprim[p]];
cline=y0-m*x0;
y=m*x+cline;
f={(x1/Cosh[Zeta0[i]])^2+(y1/Sinh[Zeta0[i]])^2==(tt*Dq[i])^2};
f1=Flatten[{x}/.Solve[f,{x}]];Print[f1];
If[f1=={},Print["No intersect!"];ntot+2,it=Times@@Dimensions[f1];If[it==1,x=f1[[1]];r1={x,y},x=f1[[1]];r1={x,y};Dis1=Abs[x+I*y-t];
Clear[x];
x=f1[[2]];r2={x,y};Dis2=Abs[x+I*y-t];
Print["Dis1EllipseIntersectPoint = ", Dis1, "    Dis2EllipseIntersectPoint = ",Dis2];
If[Dis1<=10^(-4)||Dis2<=10^(-4),If[Dis1<=10^(-4),point=r1[[1]]+I*r1[[2]],point=r2[[1]]+I*r2[[2]]],If[Dis1>Dis2,point=r2[[1]]+I*r2[[2]],point=r1[[1]]+I*r1[[2]]]]
]];
point]

Testinterface1[j_]:=Block[{i,ii,r,t,t1,t2,Dis1,Dis2,Dis3,f1,f2,f3,f4,f44,f5,sign,ltest0},  ltest0 = 0.4*(2.0*lnew + delta);
ii=IntersectCheckQ[j];
Print["IntersectCheckQ  ", ii];
i=If[Length[ii]==1,ii[[1]],ntot+2];
f1=If[i==ntot+2,Print["Crack penetration!"];{3,ii[[1]]},If[ i!=0,{r=EllipseIntersectPoint[i,ntot,tt[i]];
t2 = StartPoint[j-1];
t1 = StartPoint[j]; t = EndPoint[j];
Dis1 = Abs[Sqrt[(Re[t] - Re[t1])^2 + (Im[t] - Im[t1])^2]];
Dis2 = Abs[Sqrt[(Re[r] - Re[t1])^2 + (Im[r] - Im[t1])^2]];
Print["Dis1Testinterface1 = ", Dis1, "    Dis2Testinterface1 = ",Dis2];
Dis3 = Dis2 - Dis1;
(*Print["Dis3 = ", Dis2 - Dis1];*)
sign=If[Abs[Dis3]<ltest0,0,-1];
If[Abs[Dis3]<=10^(-9),Print["Dis3 = ", Dis3];Print["Continue Func***"];{1,i}(*;AddCrackInterface[j,i,GammaInt,GammaFrac]*),{Print["Dis3 = ", Dis3];
If[Abs[Dis3] < ltest0, f2 = 0.5 (r + t1);
  f3 = Abs[0.5 (r - t1)]; f4 = tetatemp[[j]]; f5 = tetatempprim[[j]]; 
  sign = 0, If[Sign[Dis3] == -1, f2 = 0.5 (r + t2);
   f3 = Abs[0.5 (r - t2)]; 
   f4 = myArcTan[Re[r - t2], Im[r - t2]]*180.0/Pi;
   f44 = f4 Degree;
   f5 = ArcTan[Sin[f44]/Cos[f44]]*180.0/Pi, Unevaluated[Sequence[]]]];
   If[Dis3 < ltest0, Print["f2=  ", f2, "   f3=  ", f3, "   tetatemp = ", f4, 
  "   tetatempprim = ", f5, "  sign =  ", sign];
ChangeSize[j, f3, f2, sign, f4, f5];Print["It needs reloading"];{2,sign}, {0}]}]},{0} ]]]


TestPenetration[j_]:=Block[{i,ii,r,t,t1,t2,Dis1,Dis2,Dis3,f1,f2,f3,f4,f5,f44,sign,ltest0},
(*Get["geometry.m"];
Needs["geometry`"];*)
ltest0 = 0.5*(2.0*lnew + 0delta);
ii=IntersectCheckQ[j];
Print["IntersectCheckQ  ", ii];
i=ii[[1]];
Print["i = ", i];
If[ i!=0,{r=EllipseIntersectPoint[i,ntot,tt[i]];
Print["r =  ",r, "   tt[i] =   ",tt[i]];
t2 = StartPoint[j-1];
t1 = StartPoint[j]; t = EndPoint[j];
Dis1 = Abs[Sqrt[(Re[t] - Re[t1])^2 + (Im[t] - Im[t1])^2]];
Dis2 = Abs[Sqrt[(Re[r] - Re[t1])^2 + (Im[r] - Im[t1])^2]];
Print["Dis1Penetrat = ", Dis1, "    Dis2Penetrat = ",Dis2];
Dis3 = Dis2 - Dis1;
sign=If[Abs[Dis3]<ltest0,0,-1];
Print["sign = ", sign];
If[Abs[Dis3]<=10^(-9),Print["Dis3Penetrat = ", Dis3];Print["Continue Func***"];{1,i},{Print["Dis3 = ", Dis3];
If[Abs[Dis3] < ltest0, f2 = 0.5 (r + t1);
  f3 = Abs[0.5 (r - t1)]; f4 = tetatemp[[j]]; f5 = tetatempprim[[j]]; 
  sign = 0, If[Sign[Dis3] == -1, f2 = 0.5 (r + t2);
   f3 = Abs[0.5 (r - t2)]; 
   f4 = myArcTan[Re[r - t2], Im[r - t2]]*180.0/Pi;
   f44 = f4 Degree;
   f5 = ArcTan[Sin[f44]/Cos[f44]]*180.0/Pi, Unevaluated[Sequence[]]]];
   If[Dis3 < ltest0, Print["f2=  ", f2, "   f3=  ", f3, "   tetatemp = ", f4, 
  "   tetatempprim = ", f5, "  sign =  ", sign];
ChangeSize[j, f3, f2, sign, f4, f5];Print["It needs reloading"];{2,sign}, {0}]}]},{0} ]]


PassInclusion[pp1_]:=Block[{rt,rt1,rt2,rt3,ii,ninitial1},
Print["------------------------------PassInclusion--------------------------------"];
Print[" PassInclusion =  ", pp1];
ii=IntersectCheckQ[pp1];
Print["ii PassInclusion =  ", ii];
rt3=1.9*lnew+delta;
rt2=EndPoint[pp1];
If[ii[[1]]!=0,
{rt1=EllipseIntersectPoint[ii[[1]],pp1,tt[ii[[1]]]];
rt4=Abs[EndPoint[ntot]-EllipseIntersectPoint[ii[[1]],ntot,tt[ii[[1]]]]];
rt=Abs[rt1-rt2];Print["rt: ", rt];Print["rt4: ", rt4];
If[rt>=rt3 || rt4>=rt3,
Print["Delete"];DeleteInc[ii[[1]]];
Unprotect[ninitial,lltest];
ninitial1 = ninitial;
Clear[ninitial];
ninitial = ninitial1 - 1;
writestart;
ninitial=readstart[[1]];
Print["ninitial after passInclusion",ninitial];
lltest=readstart[[2]],Print["Not passed!"]]},Print["Not passed1!"]];
]


End[] 
Protect@@Names["EllipseInterfaceCheck`*"]
EndPackage[] 

BeginPackage["EllipseInterfaceCheck`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`"}] 
EndPackage[] 













































































































