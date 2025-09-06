(* ::Package:: *)

BeginPackage["InterfaceFunc`"]
Unprotect@@Names["InterfaceFunc`*"];
ClearAll@@Names["InterfaceFunc`*"];


Energyetaeta::usage="Energyetaeta[p_] ";
geta0::usage="";
getaeta::usage="getaeta[eta_, p_]";
NewInc::usage="NewInc[p_, ThetaLinetemp_, GammaInt_, GammaFrac_]";
yprime::usage="yprime[i_, x_, y_]";
Endpointellipse::usage="Endpointellipse[i_, p_]";
AddInterfaceEllipseCrack::usage="AddInterfaceEllipseCrack[p_,i_]";
AddInterfaceLineCrack::usage="AddInterfaceLineCrack[p_,i_, ThetaNew_]";
AddCrackInterface::usage="AddCrackInterface[p_,i_,GammaInt_,GammaFrac_]";
AddCrackInterface1::usage="[i_,GammaInt_,GammaFrac_,GammaInc_]";
AddGBCrack::usage="AddGBCrack[i_,p_]";
CrackPenetrate::usage="CrackPenetrate[GammaGB_,GammaFrac_]";
CrackPenetrateGB::usage="CrackPenetrateGB[GammaGB_,GammaFrac_,i_]";
AddPenetratInc::usage="AddPenetratInc[p_,i_,GammaGB_,GammaFrac_]: New emanated crack will be added to the inclusion p on the GB i with GB fracture energy (GammaGB) and 
material fracture energy (GammaFrac).";
AddCrackpenetration::usage="AddCrackpenetration[i_,p_,emenatingPoint_]";
sigmmaLinemethod::usage = " sigmmaLinemethod[p_]: ";
sigmmaLinemethodLength::usage = " sigmmaLinemethodLength[p_,ldelta_]";
SIFLinemethodLengtht::usage = " SIFLinemethodLength[p_,ldelta_], it is useful for elleptic inclusions.";
sigmmaPointmethod::usage = " sigmmaPointmethod[p_]: ";
SigmmaC::usage = " SigmmaC[p]: It calculates the critical stress.";
ChangeLNEW2::usage = "ChangeLNEW2[p_,q_,deltat_]: It finds the appropriate length (lnew) for the (q)th inclusion with the deltet distance between pth and (q)th inclusion. This function works well just for ttest=2.";
StoredStressInct::usage="StoredStressInc[j_,ttest_]";
StoredStressInc::usage="StoredStressInc[j_,ttest_]";



Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Needs["eta0p`"];
Needs["Changesize`"];
Needs["EllipseInterfaceCheck`"];
Needs["Addone`"];
Get["Mergecrack.m"]
Needs["Mergecrack`"];



ChangeLNEW2[p_,q_,deltat_]:=Block[{etatemp,Xi1,Ztemp,lt},
etatemp=eta0[p];
Xi1=Zeta0[p]+I*etatemp;
Ztemp=dd[p]Cosh[Xi1];
lt=Abs[StartPoint[q]-Ztemp]/2.;
ChangeLNEW1[deltat,lt]]


Energyetaeta[p_] := 
 Block[{rr, rr1}, Clear[Xi, z, x, y]; rr = SigmaEtaEta[p]^2;(*+0Re[
  Tau0[1]+TauFar[1]]^2+0Im[Tau0[1]+TauFar[1]]^2*);
  Xi = XiQ[p];
  z = x + I*y;
  rr]


geta0[p_] := 
 Block[{f, rr,rsign}, rr = SigmaEtaEta[p]; Xi = Zeta0[p] + I*eta0[p]; Print["SigmaEtaEta0 =     ",rr, "  Xi =  ", Xi];rsign=Sign[rr];
  Re[rsign*rr*Conjugate[rr]]]


(*getaeta[eta_, p_] := Block[{rr, r0, r1, r2, xp, x, y, g1, r},
  rr = Energyetaeta[p];
  r0 = ((Dq[p] + r*Cos[eta])/l1[p])^2 + ((r*Sin[eta])/l2[p])^2;
  r1 = r /. Solve[r0 == 1, {r}];
  r2 = If[r1[[1]] > 0, r1[[1]], r1[[2]]];
  Print["r2 =  ", r2];
  xp = Dq[p] + r2*Cos[eta] + I*(r2*Sin[eta]);
  x = Re[xp*Exp[I*Thetaprim[p]] + zcenter[p]];
  y = Im[xp*Exp[I*Thetaprim[p]] + zcenter[p]];
  g1 = Re[rr]]*)

getaeta[eta_, p_] := Block[{rr, r0, r1, r2, xp, x, y, r,x1,y1,rsign},

rr = Energyetaeta[p];

x1 = (x - Re[zcenter[p]])*Cos[Theta[p]] + (y - Im[zcenter[p]])*
     Sin[Theta[p]];
  y1 = (y - Im[zcenter[p]])*Cos[Theta[p]] - (x - Re[zcenter[p]])*
     Sin[Theta[p]];

  r0 = ((Dq[p] + r*Cos[eta])/l1[p])^2 + ((r*Sin[eta])/l2[p])^2;
  r1 = r /. Solve[r0 == 1, {r}];
  r2 = If[r1[[1]] > 0, r1[[1]], r1[[2]]];
  Print["r2 =  ", r2];
  xp = Dq[p] + r2*Cos[eta] + I*(r2*Sin[eta]);
  x = Re[xp*Exp[I*Thetaprim[p]] + zcenter[p]];
  y = Im[xp*Exp[I*Thetaprim[p]] + zcenter[p]];
  
(*Print["z coordinate =  ",ArcCosh[(x1+I*y1)/Dq[p]]];*)

rr = SigmaEtaEta[p];
rsign=Sign[rr];
 Xi = ArcCosh[(x1+I*y1)/Dq[p]]; 
Print["SigmaEtaEta =   ",rr];
{Re[ rsign* rr*Conjugate[rr]],r2}]


(*I am not sure if I have to write ThetaNew in the Block variable list!*)
NewInc[p_, ThetaLinetemp_, GammaInt_, GammaFrac_] := 
 Block[{t,f,ff, etaa, ThetaLine, gGBdir,gGBdir1,gGBdir2, eta,eta1, getaMAX, gtemp1, gtemp2,gtemp3,gtempGB,ThetaNew,gGBdir11,gGBdir21},
 t=If[Abs[Cos[ThetaLinetemp Degree]]<=10^(-8),ThetaLinetemp,ArcTan[Sin[ThetaLinetemp Degree]/Cos[ThetaLinetemp Degree]]*180/Pi ];
 f = -Thetaprim[p]*180.0/Pi + t;
  ff = -Sign[f] 180 + f;
  Clear[eta, x, y];
  eta = f*Pi/180.0;
  Print["eta = ", eta*180.0/Pi];
  gGBdir11 = getaeta[eta, p];
  gGBdir1  = If[gGBdir11[[2]]>l1[p],-10^(9.),gGBdir11[[1]]];
  
  Clear[x, y];
  eta1 = ff*Pi/180.0;
  Print["eta = ", eta1*180.0/Pi];
  gGBdir21 = getaeta[eta1, p];
  gGBdir2  = If[gGBdir21[[2]]>l1[p],-10^(9.),gGBdir21[[1]]];
  
  
  getaMAX = geta0[p];
  gtemp1 = gGBdir1/GammaInt;
  gtemp2 = gGBdir2/GammaInt;
  gtemp3 = getaMAX/GammaFrac;
  Print["gGBdir1/GammaInt = ", gtemp1];
  Print["gGBdir2/GammaInt = ", gtemp2];
  Print["geta0/GammaFrac = ", gtemp3];
  
  Print["The minus sign in the energy means that the stress field is the minus."];
 If[Sign[gtemp1]==-1 && Sign[gtemp2]==-1 && Sign[gtemp3]==-1,Print[Style["Stress fields are negative and crack is arrested! ", FontColor -> Red]]];
  gtempGB = 
   If[gtemp1 > gtemp2, 
    ThetaNew = (eta + Thetaprim[p])*180.0/Pi; gtemp1, 
    ThetaNew = (eta1 + Thetaprim[p])*180.0/Pi; gtemp2];
  If[gtempGB > gtemp3 ,If[gtempGB<10^(-5),Print[Style["The result is surlly rong!", FontColor -> Red]]]; Print["Define new inclusion  ", ThetaNew]; 
   {1,ThetaNew,gtempGB}, If[gtemp3<10^(-5),Print[Style["The result is surlly rong! and crack is arrested.", FontColor -> Red]]];Print["AddInclusion[1]"];{0,0,gtemp3}]]


yprime[i_, x_, y_] := 
 Block[{f, f1, x1, y1, t}, 
  x1 = (x - Re[zcenter[i]])*Cos[Theta[i]] + (y - Im[zcenter[i]])*
     Sin[Theta[i]];
  y1 = (y - Im[zcenter[i]])*Cos[Theta[i]] - (x - Re[zcenter[i]])*
     Sin[Theta[i]];
  (*t = (Sin[Theta[i]] x1 + (Coth[Zeta0[i]]^2) Cos[Theta[i]] y1);
  f = If[t != 0, 
    Chop[(-Cos[Theta[i]] x1 + (Coth[Zeta0[i]]^2) Sin[
          Theta[i]] y1)/(Sin[Theta[i]] x1 + (Coth[Zeta0[i]]^2) Cos[
          Theta[i]] y1)]];
  f1 = If[t == 0, 90, ArcTan[f]*180.0/Pi]*)
f1=myArcTan[(Sin[Theta[i]] x1 + (Coth[Zeta0[i]]^2) Cos[ Theta[i]] y1),(-Cos[Theta[i]] x1 + (Coth[Zeta0[i]]^2) Sin[ Theta[i]] y1)]*180/Pi]



(*This function finds the endpoint of the pth inclusions on the interface line of the ith inclusion.*)
Endpointellipse1[i_, p_] := 
 Block[{x, y, x1, y1, t, r1, r2, r3, r4,r5,r6,v1,v2,f, tt1,j}, 
  tt1 = tt[i](*1 + delta/l1[i]*);
  x1 = (x - Re[zcenter[i]])*Cos[Theta[i]] + (y - Im[zcenter[i]])*
     Sin[Theta[i]];
  y1 = (y - Im[zcenter[i]])*Cos[Theta[i]] - (x - Re[zcenter[i]])*
     Sin[Theta[i]];
  t = EndPoint[p];
  r1 = {(x1/Cosh[Zeta0[i]])^2 + (y1/Sinh[Zeta0[i]])^2 == (tt1*
        Dq[i])^2, 
    (x - Re[t])^2.0 + (y - Im[t])^2.0 == (delta + 
       2.0*lnew)^2.(*lnew1*)};
  r2 = {x, y} /. NSolve[r1, {x, y}];
  r3 = Select[r2, Element[#, Reals] &];
 If[r3=={},Print[Style["The lnew size is larger than objects in the system!", FontColor -> Red]]];
  (*Print["r3 = ",r3];*)
 r4={Cos[Thetaprim[p ]],Sin[Thetaprim[p ]]};
r5 = Table[
    myArcTan[(r3[[j, 1]] - Re[t]), (r3[[j, 2]] - Im[t])], {j, 2}];
r6=Table[{Cos[j],Sin[j]},{j,{r5[[1]],r5[[2]]}}];
(*Print["r6   ",Re[r6]];
Print["r4   ",Re[r4]];*)
v1=Re[VectorAngle[r4,r6[[1]]]];
v2=Re[VectorAngle[r4,r6[[2]]]];
(*Print["v1   ",v1, "   v2   ", v2];*)
If[v1<v2,j = 1, j = 2];
f = r3[[j, 1]] + I*r3[[j, 2]]
  ]


(*This function finds the endpoint of the pth inclusions on the interface line of the ith inclusion.*)
Endpointellipse[i_, p_]:=Block[{tt1,t,t1,t2,t3,t4,t5,f,lines,test,p1,p2},
tt1 = tt[i];
t=EndPoint[p];
t1=Endpointellipse1[i,p];
f={{Re[t],Im[t]},{Re[t1],Im[t1]}};

t2 = EndPoint[i];
  t3 = StartPoint[i];


t4=(zcenter[i]+ l1[i]*tt1)Exp[I*Thetaprim[i]];
t5=(zcenter[i]- l1[i]*tt1)Exp[I*Thetaprim[i]];
t6=(Min[Abs[t-t4],Abs[t-t5]]-delta)/2.;


p1=f;
p2={{Re[t2],Im[t2]},{Re[t3],Im[t3]}};
lines={p1,p2};

test=segsegintersection[lines];
If[test[[1]]==True,Changelnewsize[t6,delta];Clear[t1];t1=Endpointellipse1[i,p],t1]]


(*add one crack along with the p=ntot and around the inclusion ith.*)
AddInterfaceEllipseCrack[p_,i_]:=Block[{t,t1,t2,r,r1,r2,r3},
Print["---------------AddInterfaceEllipseCrack-----------------"];
t=EndPoint[p];
t1=Endpointellipse[i,p];
r=(Re[t1]-Re[t]);
r1=myArcTan[r,(Im[t1]-Im[t])]*180.0/Pi;
r2= r1 Degree;
r3 =ArcTan[Sin[r2]/Cos[r2]]*180.0/Pi;
t2=t+lnew1*Exp[I*(r2)];
Print["theta=  ", r3,"   thetatemp=   ", r1];
ChangeSize[p, lnew, t2,+1,r3,r1]]

AddCrackInterface[p_,i_,GammaInt_,GammaFrac_]:=Block[{t,t2,r1,r2},
t=EndPoint[p];
t2=yprime[i,Re[t],Im[t]];
r1=NewInc[p,t2,GammaInt,GammaFrac];
r2 =r1[[1]];
If[r2==1,AddInterfaceEllipseCrack[p,i];{0},Print["Here I have to add AddInclusion!"];AddOne[p];{1}]
]


AddInterfaceLineCrack[p_, ThetaNew_]:=Block[{t,t2,r1,ThetaNew1,rt},
t=EndPoint[p];
t2=t+lnew1*Exp[I*( ThetaNew Degree)];
rt= ThetaNew Degree;
ThetaNew1= If[Abs[Cos[rt]]<=10^(-8), ThetaNew,ArcTan[Sin[rt]/Cos[rt]]*180.0/Pi];
ChangeSize[p, lnew, t2,+1,ThetaNew1 ,ThetaNew ]]



AddGBCrack[i_,p_,GammaInt_,GammaFrac_]:=Block[{f,xp,yp,xq,yq,ThetaNew,t,r1,r2,r3},
Get["geometry.m"];
Needs["geometry`"];
f=GBLine[i];
xp=f[[1]];yp=f[[2]];xq=f[[3]];yq=f[[4]];ThetaNew=f[[6]]*180/Pi;
t=EndPoint[p];
     r1=NewInc[p,ThetaNew,GammaInt,GammaFrac];
     r2 =r1[[1]];
     r3=If[Length[r1]!=1,r1[[2]]];
(*Print["r1  ",r1];*)
If[r2==1,AddInterfaceLineCrack[p, r3];{0},Print["Here I have to add AddInclusion!"];AddOne[p];{1}]]


AddCrackpenetration[i_,p_,emenatingPoint_]:=Block[{f,xp,yp,xq,yq,ThetaNew,t,r1,r2,r3,r4,r5,x1,y1,f11,ThetaNewt,t2,rt,ThetaNew1,tt2,x,y},
x1=(x-Re[zcenter[p]])*Cos[Theta[p]]+(y-Im[zcenter[p]])*Sin[Theta[p]];
y1=(y-Im[zcenter[p]])*Cos[Theta[p]]-(x-Re[zcenter[p]])*Sin[Theta[p]];
f11=((x1)/Cosh[Zeta0[p]])^2+((y1)/Sinh[Zeta0[p]])^2-Re[Dq[p]]^2.0;
Get["geometry.m"];
Needs["geometry`"];
f=GBLine[i];
xp=f[[1]];yp=f[[2]];xq=f[[3]];yq=f[[4]];
ThetaNew=f[[5]]*180.0/Pi;
ThetaNewt=f[[6]]*180.0/Pi;
t=emenatingPoint;
Print["emenatingPoint =  ",emenatingPoint];
r4=sigmmaLinemethod[p];
r5=SigmmaC[p];
Print["sigmmaLinemethod =    ",r4, ",   SigmmaC =   ",r5];
If[r4>=r5,Print[Style["A crack can emanate from inclusion (AddCrackpenetration[i_,p_,emenatingPoint_])!", FontColor -> Red]],Print[Style["A crack can emanate from inclusion, simulation contains error!", FontColor -> Red]]];
t2=t+lnew1*Exp[I*(ThetaNew Degree)];
rt=ThetaNew Degree;
ThetaNew1=If[Abs[Cos[rt]]<=10^(-8),ThetaNew,ArcTan[Sin[rt]/Cos[rt]]*180.0/Pi];
Clear[x,y];
Print["ThetaNew = ",ThetaNew,"  ThetaNew1 =  ", ThetaNew1,"  ThetaNewt =  ",ThetaNewt];
Print["I have to check this delta in the AddCrackpenetration function."];
(*ChangeLNEW2[p,ntot-ttest+1,delta];Get["geometry.m"];*)
tt2=t+(delta+lnew1)*Exp[I*(ThetaNewt Degree)];
x=Re[t2];y=Im[t2];
If[f11>=0,Print["f11 =  ", f11];ChangeSize[p,lnew,t2,+1,ThetaNew,ThetaNew1],Print["f11 =  ", f11];ChangeSize[p,lnew,tt2,+1,ThetaNew,ThetaNewt]]]


(*p is the inclusion that we have penetrate in, i is the GB*)
(*AddPenetratInc[p_,i_,GammaGB_,GammaFrac_]:=Block[{j,x,y,x1,t,y1,rp,Renergy,getaMAX,gtemp1,gtemp2,gtemp3,gtempGB,point,xp,yp,rr,rsign,Xi},
Clear[rp,x,y];
x1=(x-Re[zcenter[p]])*Cos[Theta[p]]+(y-Im[zcenter[p]])*Sin[Theta[p]];
y1=(y-Im[zcenter[p]])*Cos[Theta[p]]-(x-Re[zcenter[p]])*Sin[Theta[p]];

t = EndPoint[ntot];

rp=EllipseGBIntersect[i,p];
Renergy=Table[xp=rp[[j]];x=Re[xp];y=Im[xp];
(*I am not sure about this.
x=Re[xp*Exp[I*Thetaprim[p]]+zcenter[p]];
y=Im[xp*Exp[I*Thetaprim[p]]+zcenter[p]];*)
Print["Abs[t-xp] = ", Abs[t-xp]];
If[Abs[t-xp]>2.*delta,rr=SigmaEtaEta[p];
rsign=Sign[rr];
Xi=ArcCosh[(x1+I*y1)/Dq[p]];
Print["SigmaEtaEta =   ",rr, "  point=   ",xp];
{Re[rsign*rr*Conjugate[rr]]/GammaGB,xp},Unevaluated[Sequence[]]],{j,Times@Dimensions[rp][[1]]}];

getaMAX=geta0[p];
Print["Renergy =  ",Renergy];
gtemp1=Renergy[[1]][[1]];
gtemp2=If[Length[Renergy]!=1,Renergy[[2]][[1]],-10.^(15)];
gtemp3=getaMAX/GammaFrac;
Print["gGBdir1/GammaInt = ",gtemp1];
Print["gGBdir2/GammaInt = ",gtemp2];
Print["geta0/GammaFrac = ",gtemp3];

gtempGB=If[gtemp1>gtemp2,point=Renergy[[1]][[2]];gtemp1,point=Renergy[[2]][[2]];gtemp2];
If[gtempGB>gtemp3,AddCrackpenetration[i,p,point],point=dd[p]Cosh[Zeta0[p]+I*eta0[p]];AddOne[p]];
point]*)

(*New emanated crack will be added to the inclusion p on the GB i with GB fracture energy (GammaGB) and 
material fracture energy (GammaFrac).   *)
AddPenetratInc[p_,i_,GammaGB_,GammaFrac_]:=Block[{ltemp,ltemp1,j,x1,t,y1,rp,Renergy,getaMAX,gtemp1,gtemp2,gtemp3,gtempGB,point,xp,yp,rr,rsign},
Clear[rp,x,y,Xi];
x1=(x-Re[zcenter[p]])*Cos[Theta[p]]+(y-Im[zcenter[p]])*Sin[Theta[p]];
y1=(y-Im[zcenter[p]])*Cos[Theta[p]]-(x-Re[zcenter[p]])*Sin[Theta[p]];

Get["geometry.m"];
Needs["geometry`"];
ltemp=l1[ntot];
ltemp1=0.lnew/1.+ltemp;(*By putting 0 here means I am not changing the system, so it does not require Reload!*)
Print["ltemp =  ", ltemp,"   ltemp1 =  ",ltemp1];
(*Increase size of the last crack to have a penetration*)
(*The amount of 1.5 can be changed, I just checked it for my problem.*)
ChangeSize[ntot,ltemp1,zcenter[ntot],0,tetatemp[[ntot]],tetatempprim[[ntot]]];

(*Reload;*)
t=EndPoint[ntot];

rp=EllipseGBIntersect[i,p];
Renergy=Table[xp=rp[[j]];x=Re[xp];y=Im[xp];
(*I am not sure about this.x=Re[xp*Exp[I*Thetaprim[p]]+zcenter[p]];
y=Im[xp*Exp[I*Thetaprim[p]]+zcenter[p]];*)Print["Abs[Endpoint[ntot]-xp] = ",Abs[t-xp]];
If[Abs[t-xp]>2.*delta,Clear[Xi];rr=SigmaEtaEta[p];
rsign=Sign[rr];
Xi=ArcCosh[(x1+I*y1)/Dq[p]];
Print["SigmaEtaEta =   ",rr,"  point=   ",xp];
{Re[rsign*rr*Conjugate[rr]]/GammaGB,xp},Unevaluated[Sequence[]]],{j,Times@Dimensions[rp][[1]]}];

getaMAX=geta0[p];
Print["Renergy =  ",Renergy];
gtemp1=Renergy[[1]][[1]];
gtemp2=If[Length[Renergy]!=1,Renergy[[2]][[1]],-10.^(15)];
gtemp3=getaMAX/GammaFrac;
Print["gGBdir1/GammaInt = ",gtemp1];
Print["gGBdir2/GammaInt = ",gtemp2];
Print["geta0/GammaFrac = ",gtemp3];

(*Resize the last crack to its initial size.*)
(*ltemp=l1[ntot];*)
ChangeSize[ntot,ltemp,zcenter[ntot],0,tetatemp[[ntot]],tetatempprim[[ntot]]];

gtempGB=If[gtemp1>gtemp2,point=Renergy[[1]][[2]];gtemp1,point=Renergy[[2]][[2]];gtemp2];
If[gtempGB>gtemp3,AddCrackpenetration[i,p,point],point=dd[p]Cosh[Zeta0[p]+I*eta0[p]];AddOne[p]];
point]




(*NOT THIS CrackPenetrate[GammaFrac_]:=Block[{r,r1,j,i,point},
i=1;(*The first GB*)
GammaGB=GGB[i];
Get["geometry.m"];
Needs["geometry`"];
Clear[r,r1];
r=Flatten[Testinterface1[ntot]];
r1=r[[1]];
j=r[[2]];
Print["j =  ", j];
If[r1!=3,Unevaluated[Sequence[]];Print["No penetration!"],
{Print["Use crack replace and find etaMAX , ...."];
DeleteInc[ntot];
(*AppendCrackEnergy[ntot,ttest];*)
(*If we use GB or not!*)
If[GB==1,point = AddPenetratInc[j,i,GammaGB,GammaFrac],point = dd[j]Cosh[Zeta0[j]+I*eta0[j]];AddOne[j]];
Get["geometry.m"];
Needs["geometry`"];
(*routput1[2,inc];*)
ReplaceInc[j, ttest];MergeCrack[ttest+1];OutputTemp[];

PackageReload["run1`", KillShadowing -> True];
PackageReload["run1`", KillShadowing -> True];
CloseKernels[subkernel]; initsub[];
PackageReload["SIF`", KillShadowing -> True];
CloseKernels[subkernel];initsub[];
AppendCrackEnergy[ntot-1,ttest+1];}, Print["crack deflect"]];
point]*)


CrackPenetrate[GammaFrac_]:=Block[{r,r1,rt,j,i,point,ninitialtemp},
i=1;(*The first GB*)
GammaGB=GGB[i];
Get["geometry.m"];
Needs["geometry`"];
Clear[r,r1];
r=Flatten[Testinterface1[ntot]];
r1=r[[1]];
j=If[Length[r1]!=1,r[[2]]];
Print["j =  ", j];
(*Print["ChangeLNEW2[j,ntot-ttest,delta]   ", ntot];Get["geometry.m"];*)
If[r1!=3,Unevaluated[Sequence[]];Print["No penetration!"],
{Print["Use crack replace and find etaMAX , ...."];
rt=Flatten[TestPenetration[ntot]];
Print["Penetration: ", rt];
If[rt[[1]]==2&&rt[[2]]==-1, 
(*AppendCrackEnergy[ntot,ttest];*)
(*If we use GB or not!*)
If[GB==1,point = AddPenetratInc[j,i,GammaGB,GammaFrac],Reload;initsub[];point = dd[j]Cosh[Zeta0[j]+I*eta0[j]];AddOne[j]];
Get["geometry.m"];
Needs["geometry`"];
(*routput1[2,inc];*)
PlaceiToj[j, ntot-1];MergeCrack[ttest+1];OutputTemp[];

(*PackageReload["run1`",KillShadowing->True];
CloseKernels[subkernel];initsub[];
PackageReload["eta0p`",KillShadowing->True];
Needs["SIF`"];
CloseKernels[subkernel];*)
Reload;initsub[];

Unprotect[ninitial,lltest];
ninitialtemp = ninitial;
Clear[ninitial];
ninitial = ninitialtemp - 1;
writestart;
ninitial=readstart[[1]];
lltest=readstart[[2]];

PropagationTest[ntot-1,2,ttest];
AppendCrackEnergy[j,ntot-1,ttest];Print[Style["I am Here 3 (CrackPenetrate)!", FontColor -> Red],"  j=  ",j,"  ntot=  ", ntot],Unevaluated[Sequence[]]];
}];
]


CrackPenetrateGB[GammaFrac_,i_]:=Block[{r,r1,r4,r5,rt,j,point,ninitialtemp},
(*i =1; The first GB*)
GammaGB=GGB[i];
Get["geometry.m"];
Needs["geometry`"];
r=Flatten[Testinterface1[ntot]];
r1=r[[1]];
j=If[Length[r1]!=1,r[[2]]];
Print["j =  ", j];
(*ChangeLNEW2[j,ntot-ttest,delta] ;Get["geometry.m"];*)
If[r1!=3,Unevaluated[Sequence[]];Print["No penetration!"],
{Print["Use crack replace and find etaMAX , ...."];
rt=Flatten[TestPenetration[ntot]];
Print["Penetration: ", rt];
If[rt[[1]]==2&&rt[[2]]==-1,
{Print["------------------------- GB & Int ------------------------------"];
Clear[rGB1,rGB2 ,rInt1,rInt2];
Get["geometry.m"];
Needs["geometry`"];
rInt=Flatten[Testinterface1[ntot]];
rInt1=If[rInt[[1]]!=0,1,rInt[[1]]];
rInt2=If[Length[rInt]!=1,rInt[[2]],0];

Get["geometry.m"];
Needs["geometry`"];

rGB=Flatten[LineIntersectCheck[i,ntot]];
rGB1=rGB[[1]];
rGB2=If[Length[rGB]!=1,rGB[[2]],0];

If[rInt1==1&&(rGB1==1||rGB1==2),{Reload; PropagationTest[ntot,2,ttest];Print[Style["ChangeLNEW2  ", FontColor -> Blue], ntot];ChangeLNEW2[j,ntot-ttest+1,delta] ;
Print["------------------------------CrackPenetrate GB --------------------------"];
Get["geometry.m"];
Needs["geometry`"];
f=GBLine[i];
ThetaNew=f[[5]]*180/Pi;
t=EndPoint[ntot];

(*GammaFrac=  G0;*)
GammaInt= GintQ[rInt2]*EQ[rInt2];
GammaInc=GQ[rInt2]*EQ[rInt2];
Print["GammaGB =  ",GammaGB,"  GammaFrac =  ",GammaFrac, "  GammaInt =  ", GammaInt, "  GammaInc =  ", GammaInc];

     (*HERE I AM USING GAMMAINC RATHER THAN GAMMAINT BECAUSE I KNOW I HAVE CRACK PENETRATION
	IN THE GB DIRECTION.*)
(*(10^8) is for just ignoring the eta0 direction, we are aware that crack can propagate in
GB direction and it will have penetration. Just here we want to check the energy 
when a crack propagates in the GB direction with GammaInc and when it is propagating in
the interface direction with GammaInt energy.
Also, if the comparison between int and GB shows that the Inc is more likely, then it 
will automatically tests the eta0 direction.*)
	 rgb1=NewInc[ntot,ThetaNew,GammaInc,(10^8)GammaFrac];
     rgb2 =rgb1[[1]];
     rgb3=If[Length[rgb1]!=1,rgb1[[3]]];

Print["---------------------------CrackPenetrate Ellipse -----------------------"];
t2=yprime[rInt2,Re[t],Im[t]];
	 
	  rinc1=NewInc[ntot,t2,GammaInt,(10^8)GammaFrac];
	  rinc2=rinc1[[2]];
	  rinc3 =rinc1[[3]];
Print["CHECK THIS PART IF I NEED TO ADD THIS!"];
r4=sigmmaLinemethod[j];
r5=SigmmaC[j];
Print["sigmmaLinemethod =    ",r4, ",   SigmmaC =   ",r5];
If[r4>=r5,Print[Style["A crack can emanate from inclusion 1!", FontColor -> Red]],Print[Style["A crack can emanate from inclusion, simulation contains error!", FontColor -> Red]]];

If[rgb3<rinc3,jj=Interface1[rInt2,1];Print["Crack interface!"],
(*I ADD AppendCrackEnergy[j,j,ttest+1] PART, CHECK IF IT IS CAUSING ANY PROBLEM!*)
{If[GB==1,Print[Style[" ntot  ", FontColor -> Red],ntot];AppendCrackEnergy[j,j,ttest];point = AddPenetratInc[j,i,GammaGB,GammaFrac],Reload;initsub[];point = dd[j]Cosh[Zeta0[j]+I*eta0[j]];AddOne[j]];
Get["geometry.m"];
Needs["geometry`"];
PlaceiToj[j, ntot-1];MergeCrack[ttest+1];OutputTemp[];
Reload;initsub[];

Copyoutput;
Get["geometry.m"];
Needs["geometry`"];
DeleteInc[ntot];
Reload;initsub[];
AppendCrackEnergy[j,ntot,ttest-1];
OutputTemp[];Get["geometry.m"];

Unprotect[ninitial,lltest];
ninitialtemp = ninitial;
Clear[ninitial];
ninitial = ninitialtemp - 1;
writestart;
ninitial=readstart[[1]];
lltest=readstart[[2]];

(*AppendCrackEnergy[j,ntot-1,ttest];
PropagationTest[ntot-1,2,ttest];*)
Print[Style["I am Here 1!", FontColor -> Red]];
}
];},{If[GB==1,point = AddPenetratInc[j,i,GammaGB,GammaFrac],Reload;initsub[];point = dd[j]Cosh[Zeta0[j]+I*eta0[j]];AddOne[j]];
Get["geometry.m"];
Needs["geometry`"];
PlaceiToj[j, ntot-1];MergeCrack[ttest+1];OutputTemp[];

Copyoutput;
Get["geometry.m"];
Needs["geometry`"];
DeleteInc[ntot];
Reload;initsub[];
AppendCrackEnergy[j,ntot,ttest-1];
OutputTemp[];Get["geometry.m"];


Unprotect[ninitial,lltest];
ninitialtemp = ninitial;
Clear[ninitial];
ninitial = ninitialtemp - 1;
writestart;
ninitial=readstart[[1]];
lltest=readstart[[2]];

Print[Style["I am Here 2!", FontColor -> Red]];
}]

},Unevaluated[Sequence[]]]


}]]


AddCrackInterface1[i_,GammaInt_,GammaFrac_,GammaInc_]:=Block[{r,r1,r2,r3,r4,f,j},f=AddCrackInterface[ntot,i,GammaInt,GammaFrac];
Print["First Add Crack Interface"];
Get["geometry.m"];
Needs["geometry`"];
Clear[r,r1];
r=Flatten[Testinterface1[ntot]];
r1=r[[1]];
r2=r[[2]];
Print["inclusiom penetration =  ", r2];
If[r1!=3,Unevaluated[Sequence[]];Print["No penetration!"],Print["In AddCrackInterface1 func I am changing AddCrackInterface by multiplying EQ"];
{DeleteInc[ntot];
AddCrackInterface[ntot,i,GammaInt*E0,GammaInc*EQ[r2]];
Print["Second Add Crack Interface"];
r4=CrackPenetrate[GammaFrac]}];
r4
]


sigmmaLinemethod[p_] := Block[{r, Xi, z1,x1,y1,ldelta,ztemp,slop,sloptemp},
  etatemp = eta0[p];
  ztemp = Dq[p] Cosh[Zeta0[p] + I*etatemp];
  sloptemp = yprime[p, Im[ztemp], Re[ztemp]];
  slop = Tan[sloptemp*Pi/180];
  y1 = Im[ztemp] + slop*(x1 - Re[ztemp]);
  r = Re[Sigma22m[p]];
  Xi = ArcCosh[z1/Dq[p]];
  z1 = x1 + I*y1;
ldelta=(2 lnew + delta)/(1 + slop);
Print["ldelta  =  ",ldelta];
  NIntegrate[
   r, {x1, Re[ztemp],ldelta  + Re[ztemp]}]/ldelta
  ]

sigmmaLinemethodLength[p_,ldelta_]:= Block[{r,r12,f22,f12,f, Xi, z1,x1,y1,ztemp,slop,sloptemp},
  etatemp = eta0[p];
  ztemp = Dq[p] Cosh[Zeta0[p] + I*etatemp];
  sloptemp = yprime[p, Im[ztemp], Re[ztemp]];
  slop = Tan[sloptemp*Pi/180];
  y1 = Im[ztemp] + slop*(x1 - Re[ztemp]);
  r = Re[Sigma22m[p]];
  r12=Re[Sigma12m[p]];
  Xi = ArcCosh[z1/Dq[p]];
  z1 = x1 + I*y1;
  Print["ldelta  =  ",ldelta];
  f22= NIntegrate[
   r, {x1, Re[ztemp],ldelta  + Re[ztemp]}]/ldelta;
  f12= NIntegrate[
   r12, {x1, Re[ztemp],ldelta  + Re[ztemp]}]/ldelta;
	f=f22+I*f12
  ]


SIFLinemethodLengtht[p_,ldelta_]:= Block[{etatemp,r,r12,f22,f12,f, Xi, z1,x1,y1,ztemp,slop,sloptemp},
  etatemp = eta0[p];
  ztemp = Dq[p] Cosh[Zeta0[p] + I*etatemp];
  sloptemp = yprime[p, Im[ztemp], Re[ztemp]];
  slop = Tan[sloptemp*Pi/180];
  y1 = Im[ztemp] + slop*(x1 - Re[ztemp]);
Clear[Xi];
  r = Re[Sigma22m[p]];
  r12=Re[Sigma12m[p]];
  Xi = ArcCosh[z1/Dq[p]];
  z1 = x1 + I*y1;
  Print["ldelta  =  ",ldelta];
  f22= NIntegrate[
   Sqrt[(x1- Re[ztemp])*2.Pi]r, {x1, Re[ztemp],ldelta  + Re[ztemp]}]/ldelta;
  f12= NIntegrate[
   Sqrt[(x1- Re[ztemp])*2.Pi]r12, {  x1, Re[ztemp],ldelta  + Re[ztemp]}]/ldelta;
	f=f22+I*f12
  ]


sigmmaPointmethod[p_] := Block[{r, Xi, z1, y1, x1},
  etatemp = eta0[p];
  ztemp = Dq[p] Cosh[Zeta0[p] + I*etatemp];
  sloptemp = yprime[p, Im[ztemp], Re[ztemp]];
  slop = Tan[sloptemp*Pi/180];
  y1 = Im[ztemp] + slop*(x1 - Re[ztemp]);
  r = Re[Sigma22m[p]];
  Xi = ArcCosh[z1/Dq[p]];
  z1 = x1 + I*y1;
  x1 = (lnew/2 + delta/4)/(1 + slop) + Re[ztemp];
  r]


(*Stored stress in the inclusion j due to the presence of crack with ttest number of micro-cracks.*)
StoredStressInc[j_,ttest1_]:=Block[{f,f1,f2,r,r1,rr0,rr1,rr2,rr},
Clear[Xi];
(*Print["I HAVE CHANGED {ii,ntot-ttest+1,ntot} TO {ii,ntot-ttest+1,ntot-1}, TO ELEMINATE THE EFFECT OF THE LAST MICROCRACK ON THE INCLUSION J."];*)
f  =  Table[If[ii!=j,ii,Unevaluated[Sequence[]]],{ii,Min[ntot-ttest1+1,ntot-1],Max[ntot-ttest1+1,ntot-1]}];
Print["coefficient test=  ", f];
f1 = Sum[means22[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];
f2 = Sum[means12[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];

r=SIFLinemethodLengtht[j,0.01l1[j]]; (*It is a distance where SIF is averaged over it.*)
rr=SigmaEtaEta[j];
Xi=Zeta0[j]+I*eta0[j];
Print["SigmaEtaEta[j] =  ",rr];
Print["These Energy and SIF are not necessarily correct."];
Print["r = ",r];
r1 = f1 + I*f2 +r]


StoredStressInct[j_,ttest1_]:=Block[{f,f1,f2,r,r1,rr0,rr1,rr2,rr},
Clear[Xi];
(*Print["I HAVE CHANGED {ii,ntot-ttest+1,ntot} TO {ii,ntot-ttest+1,ntot-1}, TO ELEMINATE THE EFFECT OF THE LAST MICROCRACK ON THE INCLUSION J."];*)
f  =  Table[If[ii!=j,ii,Unevaluated[Sequence[]]],{ii,Min[ntot-ttest1+1,ntot-1],Max[ntot-ttest1+1,ntot-1]}];
Print["coefficient test=  ", f];
f1 = Sum[means22[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];
f2 = Sum[means12[j, jj]*Sqrt[Pi*l1[j]],{jj,f}];

rr0 = ((Dq[j] + r*Cos[eta0[j]])/l1[j])^2 + ((r*Sin[eta0[j]])/l2[j])^2;
  rr1 = r /. Solve[rr0 == 1, {r}];
  rr2 = If[rr1[[1]] > 0, rr1[[1]], rr1[[2]]];
  Print["r =  ", rr2];

eta0temp=eta0[j];
rs=sigmmaLinemethodLength[j,2.l1[j]e[j]^2];
r=Sqrt[2*Pi*(dd[j](Cosh[Zeta0[j]+I*eta0temp]-1))]rs;
rr=SigmaEtaEta[j];
Xi=Zeta0[j]+I*eta0temp;
Print["SigmaEtaEta[j] =  ",rr];
Print["rs   ", rs];
Print["r = ",r];
Print["rCoeff = ",Sqrt[2*Pi*(dd[j](Cosh[Zeta0[j]+I*eta0temp]-1))]];
r1 = f1 + I*f2 +r]


EModulus=2.*Mu0(1+Nu0)
SigmmaC[i_]:=Sqrt[EModulus*G0/2./l1[i]/1000](*GPa*)


End[] 
Protect@@Names["InterfaceFunc`*"]
EndPackage[] 

BeginPackage["InterfaceFunc`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","Addone`","EllipseInterfaceCheck`","Mergecrack`","EnergyTest`","LineInterfaceCheck`","CombinDCDTEMP`"}] 
EndPackage[] 






















































































































