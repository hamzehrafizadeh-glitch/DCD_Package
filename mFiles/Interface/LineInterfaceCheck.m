(* ::Package:: *)

BeginPackage["LineInterfaceCheck`"]
Unprotect@@Names["LineInterfaceCheck`*"];
ClearAll@@Names["LineInterfaceCheck`*"];



LineIntersectPoint::usage="LineIntersectPoint[i_, xp_, yp_, xq_, yq_] ";
GBLine::usage="GBLine[i_]: This gives i-th GBLine information as {xp, yp, xq, yq, ThetaLine}";
LineIntersectCheck::usage="LineIntersectCheck[i_,p_]: check intersect of the GB[i] and crack p";
EllipseGBIntersect::usage="EllipseGBIntersect[p_, i_]: This function finds the intersection points of the p GB by the i-th \
inclusion. file LineInterfaceCheck";
segsegintersection::usage="";


Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Needs["eta0p`"];
Needs["Changesize`"];
Needs["InterfaceFunc`"]
Needs["EllipseInterfaceCheck`"]



GBLine[i_]:=Block[{f,f1,f2,f3,f4,f5,xp,xq,yp,yq},
f=Flatten[vertexList[[i]]];
xp=f[[1]];
yp=f[[2]];
xq=f[[3]];
yq=f[[4]];
f1=ArcTan[xq - xp,yq - yp];
f2=If[Abs[Cos[f1]]<=10^(-8),f1,ArcTan[Sin[f1]/Cos[f1]]];
f3=myArcTan[-xq + xp,-yq + yp];
f4=myArcTan[xq - xp,yq - yp];
f5=If[f3==f2,f4,f3];
Flatten[{f,f2,f5}]]


(*LineIntersectPoint[i_, xp_, yp_, xq_, yq_] := 
 Block[{x, y, r, x1, y1, t, x0, y0, m, cline, f, f1,f2,f3, ltest0,t1,t2,Dis1,Dis2,Dis3,sign, rtest,ThetaNew,ff1}, 
  x1 = (x - Re[zcenter[i]])*Cos[Theta[i]] + (y - Im[zcenter[i]])*
     Sin[Theta[i]];
  y1 = (y - Im[zcenter[i]])*Cos[Theta[i]] - (x - Re[zcenter[i]])*
     Sin[Theta[i]];
  t = EndPoint[i];
  t1 = StartPoint[i];
  t2 = StartPoint[i-1];
  (*ff1=ArcTan[xq - xp,yq - yp];
  ThetaNew=If[Abs[Cos[ff1]]<=10^(-8),ff1,ArcTan[Sin[ff1]/Cos[ff1]]];*)
  x0 = Re[t]; y0 = Im[t];
  m = Tan[Thetaprim[i]];
  cline = y0 - m*x0;
  y = m*x + cline;
  f = If[xq == 
     xp, {x == xp}, {y - yp == ((yq - yp)/(xq - xp)) (x - xp)}];Print["f  ", f];
  f1 = Flatten[{x} /. Solve[f, {x}]];
   If[f1 == {}, Print["No intersect!"]; ntot + 2, x = f1[[1]];
   r = {x, y};
   (*Print["r  ",r];*)
   Dis1 = Sqrt[(Re[t] - Re[t1])^2 + (Im[t] - Im[t1])^2];
   Dis2 = Sqrt[(r[[1]] - Re[t1])^2 + (r[[2]] - Im[t1])^2];
   Dis3 = Dis2 - Dis1;
Print["Dist1=   ",Dis1,"   Dis2 =   ",Dis2];
   sign=If[Abs[Dis3]<ltest0,0,-1];
   ltest0 = 0.5*(2.0*lnew + delta);
   rtest=ThetaNew+myArcTan[-(r[[1]] - Re[t1]) ,-(r[[2]] - Im[t1])]*180/Pi ;
   (*If[ArcSin[Sin[rtest]]>10^(-8),Print["No Intersect!"];{0},*)
If[Abs[Dis3]<=10^(-9),Print["Dis3 = ", Dis3];Print["Continue Func***"];{1,i}(*here i should be changed to what ever we name as 
   the GB information*),
Print["Dis3 = ", Dis3]; 
(* I have to use this one if I work out for two directions.
 If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0,If[Sign[Dis3]==-1,f2=0.5 (r[[1]] + I*r[[2]] + t2);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t2)],Unevaluated[Sequence[]]]];*)

If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0,Unevaluated[Sequence[]]];

If[Abs[Dis3] < ltest0, Print["f2=  ",f2,"   f3=  ", f3,"  sign =  ",sign];
ChangeSize[i, f3, f2,sign,tetatemp[[i]],tetatempprim[[i]]];Print["It needs reloading"];{2,sign}, {0}]
]],{0} ]*)
(*1 = on the GB
  0 = No GB
  2 = resized and needs reloading. *)


LineIntersectPointtest[i_, xp_, yp_, xq_, yq_] := 
 Block[{ltest0,ti,x1,y1,t,t1,t2,t3,x0,y0,m,cline,f11,f12,ff,ff1,ff2,ff3,endpointtest1,endpointtest,startpointtest1,startpointtest,linetest,test,f1,f2,f3,r,Dis1,Dis2,Dis3,sign},

ltest0 = 0.5*(2.0*lnew + delta);ti=1.0+ltest0/l1[i];
x1 = (x - Re[zcenter[i]])*Cos[Theta[i]] + (y - Im[zcenter[i]])*Sin[Theta[i]];
y1 = (y - Im[zcenter[i]])*Cos[Theta[i]] - (x - Re[zcenter[i]])*Sin[Theta[i]];

 t = EndPoint[i];
 t1 = StartPoint[i];
 t2 = StartPoint[i-1];
 t3=EndPoint[i-1];
  x0 = Re[t]; y0 = Im[t];
  m = Tan[Thetaprim[i]];
  cline = y0 - m*x0;

f11={((x1)/Cosh[Zeta0[i]])^2+((y1)/Sinh[Zeta0[i]])^2==Re[ti*Dq[i]]^2.0};

f12 = If[xq == 
     xp, {x == xp}, {y - yp == ((yq - yp)/(xq - xp)) (x - xp)}];
ff=Join[f11,f12];

ff1=Solve[ff,{x,y}];
ff2={x,y}/.ff1;
ff3=Select[ff2,Element[#,Reals]&];

     Clear[x,y];
endpointtest1=f12;
x=Re[t];y=Im[t];
endpointtest=If[endpointtest1=={True},1,0];
Clear[x,y];
startpointtest1=f12;
x=Re[t3];y=Im[t3];
startpointtest=If[startpointtest1=={True},1,0];
linetest=If[endpointtest==0 &&startpointtest==1,1,0];
Clear[x,y];

test=If[ff3=={} || linetest==1 ,Print["No Intersect"];{0},
y = m*x + cline;

  f1 = Flatten[{x} /. Solve[f12, {x}]];Print["f1  ",f1];
   If[f1 == {}, Print["No intersect!"]; ntot + 2, x = f1[[1]];
   r = {x, y};
   (*Print["r  ",r];*)
   Dis1 = Sqrt[(Re[t] - Re[t1])^2 + (Im[t] - Im[t1])^2];
   Dis2 = Sqrt[(r[[1]] - Re[t1])^2 + (r[[2]] - Im[t1])^2];
   Dis3 = Dis2 - Dis1;
Print["Dist1=   ",Dis1,"   Dis2 =   ",Dis2];
   sign=If[Abs[Dis3]<ltest0,0,-1];
   (*ltest0 = 0.5*(2.0*lnew + delta);*)
   If[Abs[Dis3]<=10^(-9),Print["Dis3 = ", Dis3];Print["Continue Func***"];{1,i}(*here i should be changed to what ever we name as 
   the GB information*),
Print["Dis3 = ", Dis3]; 
 If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0,If[Sign[Dis3]==-1,f2=0.5 (r[[1]] + I*r[[2]] + t2);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t2)],Unevaluated[Sequence[]]]];

(*If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0,Unevaluated[Sequence[]]];*)

If[Dis3 < ltest0, Print["f2=  ",f2,"   f3=  ", f3,"  sign =  ",sign];
ChangeSize[i, f3, f2,sign,tetatemp[[i]],tetatempprim[[i]]];Print["It needs reloading"];{2,sign}, {0}]
],{0} ]]]


(*1 = on the GB
  0 = No GB
  2 = resized and needs reloading. *)


LineIntersectChecktest[i_,p_]:=Block[{f,xp,yp,xq,yq,ThetaNew},
f=GBLine[i];
xp=f[[1]];yp=f[[2]];xq=f[[3]];yq=f[[4]];ThetaNew=f[[5]]*180.0/Pi;
(*Print["ThetaNew=   ",ThetaNew];*)

If[Abs[(Theta[p]*180/Pi-ThetaNew )]<10^(-8)||Abs[(Theta[p]*180/Pi-ThetaNew+180)]<10^(-8),Print["Parallel and Continue Func***"];{1},Flatten[LineIntersectPoint[p, xp, yp, xq, yq]]]
]

(*This function finds the intersection points of the p GB by the i-th \
inclusion.*)




segsegintersection[lines_]:=Module[{md=Subtract@@(Plus@@#&/@lines),sub=Subtract@@#&/@lines,det},det=-Det[sub];Print["segsegintersection Det",det];
If[Abs[det]<=10^(-5),Print["Parallel and I have to check distance"];{False,0},
If[And@@(Abs[#]<=1&/@#),{True,(Plus@@#[[1]]-Subtract@@#[[1]] Last@#[[2]])/2&@{First@lines,#}},{False,1}]&@(Det[{#[[1]],md}]/det&/@({#,Reverse@#}&@sub))]];


distance[{start_,end_},pt_]:=Module[{param=((pt-start) . (end-start))/Norm[end-start]^2},EuclideanDistance[pt,start+Clip[param,{0,1}] (end-start)]];
(*Min[distance[#,p2[[1]]]&/@{p1}]
  Min[distance[#,point]&/@lines]*)


LineIntersectPoint[p_,vertexList_,ltest0_] := 
 Block[{p1,p2,lines,test,test1,r,f2,f3,t0,t,t1,t2,t3,Dis1,Dis2,Dis3,sign,test2,i},
i=1;

(*ltest0 = 0.5*(2.0*lnew + delta);*)
t0=zcenter1[[p]]+(l1[p]+ltest0)*Exp[I*Thetaprim[p]];
t = EndPoint[p];
  t1 = StartPoint[p];
 t2 = StartPoint[p-1];
 t3=EndPoint[p-1];

p1=vertexList;
p2={{Re[t1],Im[t1]},{Re[t0],Im[t0]}};
lines={p1,p2};
test=segsegintersection[lines];
Print["test =  ", test];
test1=If[test[[1]]==False,test2=Min[distance[#,p2[[1]]]&/@{p1}];Print["distance = ", test2];If[test2<=10^(-2),Print["Parallel"];{1},Print["No Intersect"];{0}],
r=test[[2]];Print["r =  ",r];
   Dis1 = Sqrt[(Re[t] - Re[t1])^2 + (Im[t] - Im[t1])^2];
   Dis2 = Sqrt[(r[[1]] - Re[t1])^2 + (r[[2]] - Im[t1])^2];
   Dis3 = Dis2 - Dis1;
Print["Dist1 LineIntersect=   ",Dis1,"   Dis2 LineIntersect =   ",Dis2];
 ltest0 = 0.5*(2.0*lnew + delta);
   sign=If[Abs[Dis3]<ltest0,0,-1];
  
   If[Abs[Dis3]<=10^(-2),Print["Dis3 LineIntersect = ", Dis3];Print["Continue Func***"];{1,i}(*here i should be changed to what ever we name as 
   the GB information*),
Print["Dis3 = ", Dis3]; 
 If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0;
    r1=tetatempprim[[p]];
	r3=tetatemp[[p]],If[Sign[Dis3]==-1,f2=0.5 (r[[1]] + I*r[[2]] + t2);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t2)];
	t4=-t2+(r[[1]] + I*r[[2]]  );
   r1=myArcTan[Re[t4],Im[t4]]*180.0/Pi;
    r2= r1 Degree;
    r3 =ArcTan[Sin[r2]/Cos[r2]]*180.0/Pi;,Unevaluated[Sequence[]]]];

(*If[Abs[Dis3]< ltest0,f2=0.5 (r[[1]] + I*r[[2]] + t1);
    f3 = Abs[0.5 (r[[1]] + I*r[[2]] - t1)];sign=0,Unevaluated[Sequence[]]];*)

If[Dis3 < ltest0, Print["f2=  ",f2,"   f3=  ", f3,"  sign =  ",sign];
ChangeSize[p, f3, f2,sign,r3,r1];Print["It needs reloading"];{2,sign}, {0}]
],{0} ]]



LineIntersectCheck[i_,p_]:=Block[{f,f1},
f=vertexList[[i]];
f1=LineIntersectPoint[p,f,0.5*(2.0*lnew + delta)]]


EllipseGBIntersect[p_, i_] := 
 Block[{f, xp, yp, xq, yq, ThetaNew, ltest0, ti, x1, y1, f11, f12, ff,
    ff1, ff2, ff3, ff4,x,y},
  Clear[x,y];
  f = GBLine[p];
  xp = f[[1]]; yp = f[[2]]; xq = f[[3]]; yq = f[[4]]; 
  ThetaNew = f[[5]]*180.0/Pi;
  
  ltest0 = 0.5*(2.0*lnew + delta); ti = 1.0 ;+ ltest0/l1[i];
  x1 = (x - Re[zcenter[i]])*Cos[Theta[i]] + (y - Im[zcenter[i]])*
     Sin[Theta[i]];
  y1 = (y - Im[zcenter[i]])*Cos[Theta[i]] - (x - Re[zcenter[i]])*
     Sin[Theta[i]];
  
  
  f11 = {((x1)/Cosh[Zeta0[i]])^2 + ((y1)/Sinh[Zeta0[i]])^2 == 
     Re[ti*Dq[i]]^2.0};
  
  f12 = If[
    xq == xp, {x == xp}, {y - yp == ((yq - yp)/(xq - xp)) (x - xp)}];
  ff = Join[f11, f12];
  
  ff1 = Solve[ff, {x, y}];
  ff2 = {x, y} /. ff1;
  ff3 = Select[ff2, Element[#, Reals] &];
  ff4 = Table[
    ff3[[j]][[1]] + I*ff3[[j]][[2]], {j, Times @@ Dimensions[ff3]/2}]
  ]


End[] 
Protect@@Names["LineInterfaceCheck`*"]
EndPackage[] 

BeginPackage["LineInterfaceCheck`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" , "matrix`" , "farfield`" ,"linearSolve`", "PhiPsi`","StressFields`","DisplacementFields`","SIF`","eta0p`","Changesize`","EllipseInterfaceCheck`"}] 
EndPackage[] 




















































