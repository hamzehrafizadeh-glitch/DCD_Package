(* ::Package:: *)

BeginPackage["Addone`"]
Unprotect@@Names["Addone`*"];
ClearAll@@Names["Addone`*"];







AddOne::usage="This function adds 1 inclusion with a specific length lnew1 to the system. Adding process is started from the last inclusion ntot-th in the input file.";
AddOneSide::usage="";
AddOneDouble::usage="";
NewCenter::usage="Center of new inclusion to add to the system as crack propagates.";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Needs["eta0p`"];




(*It obtaines the position of the new coordinate system origin where the new inclusion is
originated from the pth inclusion.*)

rtest[p_,feta_]:=Block[{f,d},Clear[x,y];d=Re[dd[p]*Exp[-I Theta[p]]];
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


NewCenter[p_,feta_]:=Block[{f1,r,r1,f,xx1,yy1,x0,y0},r=rtest[p,feta];
x0=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[feta];
y0=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[feta];

If[Dimensions[r]=={1,2},xx1=r[[1]][[1]];yy1=r[[1]][[2]],Table[s=If[Sign[r[[i]][[1]]]==Sign[x0]&&Sign[r[[i]][[2]]]==Sign[y0],1,0];
If[s==1,xx1=r[[i]][[1]];yy1=r[[i]][[2]]],{i,2}]];
f=xx1+I*yy1]
(*NewCenter[ntot] is equal to (xx1+I*yy1)*)


AddOneMaxTang[p_]:=Block[{r,x0,y0,xx1,yy1,rtt,rt,rt1,feta,rt0},feta=eta0[p];
Get["geometry.m"];
Needs["geometry`"];
x0=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[feta];
y0=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[feta];

r=NewCenter[p,feta];
(*r=NewCenter[p,0];*)
xx1=Re[r];
yy1=Im[r];

rt0=Re[N[myArcTan[(xx1-x0),(yy1-y0)]/\[Degree]]];
Print["rt0  ", rt0, "  eta0[p]=  ", feta];
rt1=If[-Pi/2<=feta<Pi/2.0,Thetaprim[p]*180/Pi+rt0,(Thetaprim[p]-Pi)*180/Pi+rt0];
Print[" rt1  ", rt1];
(*rt=If[0.00001<=rt1<=270.0 || -90.0<=rt1<=10^(-3.0) ,rt1,rt1-360];*)
rt = If[0.00001 <= rt1 <= 270.0 || -90.0 <= rt1 <= 10^(-3.0), rt1, 
  If[Abs[rt1] >= 360, 360 - Mod[rt1, 360], 
   If[Sign[rt1] == -1, rt1 + 360, rt1 - 360]]];

rtt= rt Degree;
Print["Here we might have a problem with angle!"];
Print["rt is an inclination angle: ", rt];
Thetanew = ArcTan[Sin[rtt]/Cos[rtt]]*180.0/Pi;
(*If[rt>90.0001,rt-180,If[rt<-90.01,180+rt,rt]];*)
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




AddOneSIF[p_] := Block[{r, x0, y0, xx1, yy1, rtt, rt, rt1, rt0},
  feta = SIF`eta0SIF[p]+Thetaprim[p];
  Thetanew = ArcTan[Sin[feta]/Cos[feta]]*180.0/Pi;
  Print["Thetanew SIF is: ", Thetanew];
  ThetaNewtemp = feta*180./Pi;
  zcenternew = EndPoint[p] + (lnew + delta)*Exp[I*feta];
 Print["ZcenterNew: ", zcenternew];
  
  
  tetatemp1 = Append[tetatemp, Thetanew];
  tetatemp1prim = Append[tetatempprim, ThetaNewtemp];
  l1temp1 = Append[l1temp, lnew];
  l2temp1 = Append[l2temp, lnew/l2ratio];
  NuQ1 = Append[NuQtemp, 0.0];
  MuQ1 = Append[MuQtemp, 0.0];
  zcentertemp = Append[zcenter1,zcenternew];
  zcenternewtemp1 = zcentertemp;
  ntot1 = ntot + 1; r3[]]


AddOne[p_] := If[num1==1,AddOneMaxTang[p],AddOneSIF[p],AddOneMaxTang[p]]


AddOneSide[p_,feta_]:=Block[{r,x0,y0,xx1,yy1,rtt,rt,rt1,rt0},
Get["geometry.m"];
Needs["geometry`"];
x0=Re[dd[p]*Exp[-I*Theta[p]]]*Cosh[Zeta0[p]]*Cos[feta];
y0=Re[dd[p]*Exp[-I*Theta[p]]]*Sinh[Zeta0[p]]*Sin[feta];

r=NewCenter[p,feta];
xx1=Re[r];
yy1=Im[r];

rt0=Re[N[myArcTan[(xx1-x0),(yy1-y0)]/\[Degree]]];
Print["rt0  ", rt0, "  eta0[p]=  ", feta];
rt1=If[-Pi/2<=feta<Pi/2.0,Thetaprim[p]*180/Pi+rt0,(Thetaprim[p]-0Pi)*180/Pi+rt0];
Print[" rt1  ", rt1];
rt=If[0.00001<=rt1<=270.0 || -90.0<=rt1<=10^(-3.0) ,rt1,rt1-360];
rtt= rt Degree;
Print["Here we might have a problem with angle!"];
Print["rt is an inclination angle: ", rt];
Thetanew = ArcTan[Sin[rtt]/Cos[rtt]]*180.0/Pi;
(*If[rt>90.0001,rt-180,If[rt<-90.01,180+rt,rt]];*)
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
    ntot1 = ntot + 1;r3[];]



AddOneR[p_]:={AddOneSide[p,eta0[p]]}
AddOneL[p_]:={AddOneSide[p,eta0L[p]]}

AddOneDouble[p_]:={ntemp=p-1;
If[ntemp==0,AddOne[p];AddOneL[p],AddOne[p-1];AddOne[p]]}


End[] 
Protect@@Names["Addone`*"]
EndPackage[] 

BeginPackage["Addone`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","CombinMMBE`","eta0p`","SIF`"}] 
EndPackage[] 








































