(* ::Package:: *)

 BeginPackage[ "GenerateMesh`"];
 Unprotect@@Names["GenerateMesh`*"];
 ClearAll@@Names["GenerateMesh`*"];


NTOTAL::usage ="Total number of nodes in the mesh.";
N1::usage ="It is the number of nodes at the up & Down width.";
N2::usage ="It is the number of nodes at the 2 heights.";
n1::usage ="It is the number of nodes at the up width.";
n2::usage ="It is the number of nodes at the heights.";
ncorner::usage ="It is the number of nodes at the up width.";


 Begin[ "Private`"];
 




(*------------    DISCRITIZE ---------------------*)
(*New scheme of \
discretization for non-sharp corner rectangular*)
L1 =120; L2 = 90.; 
(*Rectangular dimention*)

e1 = 1*0.06/2;
e2 = 1*0.06/2;
a3 =1;
dep = -2*e2 - L2;

Offset1 = -0;
Offset2 = 0;
n1 = (L1/a3 + 1);
n2 = (L2/a3 + 1);

N1 = 2*n1; (*It is the number of nodes at the up & Down width*)

N2 = 2*n2; (*It is the number of nodes at the 2 heights*)

ntotal = N1 + N2;
(*Here I define number of elements at each corner with respect to a3.
The idea is to have a less mismatch between elements seperation.*)
For[k = 1, k < ntotal + 1, k++,
 ncorner = k; If[\[Pi]/2/(ncorner)*e2/a3 <= 1, Break[]]]


(*This means I have changed it to not having corner correction.*)
ncorner=0;


NTOTAL = ntotal + ncorner*4;
XY = Table[{0, 0, 0, 0, 0, 0, 0}, {NTOTAL}];
Do[XY[[i]] = {Offset1 + e1 - a3 + a3*i, 0 - Offset2, a3, 0, a3, 
    Offset1 + e1 - a3 + a3*i, -Offset2}, {i, 1, N1/2}];

Do[XY[[i]] = {L1 + 0*e1 + Offset1, -Offset2 + a3 - e2 - a3*(i - n1), 
    a3, -\[Pi]/2, a3, 
    L1 + 0*e1 + Offset1, -Offset2 + a3 - e2 - a3*(i - n1)}, {i, 
   n1 + 1, n1 +n2}];

Do[XY[[i]] = {Offset1 + 2e1 - a3 +L1- (i - n1-n2-1)*a3, dep - Offset2+2e1, 
    a3, \[Pi], a3, Offset1 + 2e1 - a3 +L1- (i - n1-n2-1)*a3, 
    dep - Offset2}, {i, n1 +n2 + 1, 2n1 +n2}];

Do[XY[[i]] = {0 + Offset1, -Offset2 + a3 - 2e2 -L1+ a3*(i -( 2n1 +n2+1)), 
    a3, \[Pi]/2, 
    a3, +Offset1, -Offset2 + a3 - 2e2 -L1+ a3*(i -( 2n1 +n2+1))}, {i,  2n1 +n2+1,N1+N2}];

O1 = {+Offset1 + L1 + e1, -Offset2 - e2};
O2 = {Offset1 + L1 + e1, dep + e2 - Offset2};
O3 = {e1 + Offset1, -L2 - e2 - Offset2};
O4 = {e1 + Offset1, -e2 - Offset2};
f1 = \[Pi]/2;
f2 = 0;
f3 = -\[Pi]/2;
f4 = -\[Pi];

ntemp=ncorner+1;
(*wfunc[j_]:=2*Sin[Pi-Pi*((ncorner-i)/ncorner)]^2;*)

Do[XY[[ntotal + i, 1]] = O1[[1]] + e1*Cos[f1 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + i, 2]] = O1[[2]] + e2*Sin[f1 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + i, 3]] = \[Pi]/2/ntemp*e1;
  XY[[ntotal + i, 4]] = \[Pi]/2 - f1 - i*(\[Pi]/2)/ntemp, {i, 1, 
  ncorner}]


Do[XY[[ntotal + ncorner + i, 1]] = 
  O2[[1]] + e1*Cos[f2 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + i + ncorner, 2]] = 
  O2[[2]] + e2*Sin[f2 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + ncorner + i, 3]] = \[Pi]/2/ntemp*e1;
  XY[[ntotal + ncorner + i, 4]] = \[Pi]/2 - f2 - 
   i*(\[Pi]/2)/ntemp, {i, 1, ncorner}]

Do[XY[[ntotal + 2*ncorner + i, 1]] = 
  O3[[1]] + e1*Cos[f3 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + i + 2*ncorner, 2]] = 
  O3[[2]] + e2*Sin[f3 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + 2*ncorner + i, 3]] = \[Pi]/2/ntemp*e1;
  XY[[ntotal + 2*ncorner + i, 4]] = \[Pi]/2 - f3 - 
   i*(\[Pi]/2)/ntemp, {i, 1, ncorner}]

Do[XY[[ntotal + 3*ncorner + i, 1]] = 
  O4[[1]] + e1*Cos[f4 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + i + 3*ncorner, 2]] = 
  O4[[2]] + e2*Sin[f4 - i*(\[Pi]/2)/ntemp];
  XY[[ntotal + 3*ncorner + i, 3]] = \[Pi]/2/ntemp*e1;
  XY[[ntotal + 3*ncorner + i, 4]] = \[Pi]/2 - f4 - 
   i*(\[Pi]/2)/ntemp, {i, 1, ncorner}]

Do[XY[[i, 6]] = XY[[i, 1]]; XY[[i, 7]] = XY[[i, 2]];
  XY[[i, 5]] = XY[[i, 3]], {i, 1, NTOTAL}]



XY >> XY.dat;
ListPlot[Table[{XY[[i, 1]], XY[[i, 2]]}, {i, 1, NTOTAL}](*, 
 Axes -> None*) ,AspectRatio->Automatic]
NTOTAL


 End[];
 Protect@@Names["GenerateMesh`*"]
 EndPackage[]



























