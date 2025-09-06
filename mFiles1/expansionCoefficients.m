(* ::Package:: *)

BeginPackage["ExpansionCoefficient`"]
Unprotect @@ Names["ExpansionCoefficient`*"];
ClearAll @@ Names["ExpansionCoefficient`*"];

EtaFarPQ::usage=" EtaFarPQ = \[Eta]far[n, m, p, q]"; 
MuFarPQ::usage=" MuFarPQ = \[Mu]far[n, m, p, q]"; 

EtaClosePQ::usage=" EtaClosePQ = \[Eta][n, m, p, q] using different expansion for close interactions."; 
MuClosePQ::usage=" MuClosePQ = \[Mu][n, m, p, q] using different expansion for close interactions."; 

EtaIntgPQ::usage="EtaIntgPQ = \[Eta][n, m, p, q] using integral format."; 
MuIntgPQ::usage="MuIntgPQ = \[Mu][n, m, p, q] using integral format."; 

(*Eta::usage=" Eta"; 
Mu::usage="Mu"; *)
zpq::usage=" zpq=relation between the local coordinate of different systems"; 
ltest::usage=" ltest"; 
nmax1::usage=" nmax1";
nmax::usage=" nmax";


Begin["`Private`"] 


(*Change the path into mFiles path.*)
AppendTo[$Path,FileNameJoin[{FilePath, "mFiles"}]];
$Path

(****Get["geometry.m"];*)
Needs["geometry`"];

nmax=10
nmax1= 2 nmax+n+1

ll1 = 2 l1[1]
ltest = Sqrt[3] dd[1]  (*CHANGE THIS PART LATER*)
ll2 = 2 l2[1]
 
(*zcenter[p_] := 
 Module[{f, i, j}, i = n1 FractionalPart[(p - 1)/n1]; 
  j = IntegerPart[(p - 1)/n1];
  f = i ll1 + I j ll2  ]*)
 
zpq[p_, q_] := Module[{f},
  f = zcenter[q] - zcenter[p]]
  

  
vpq[p_, q_] := 
 Module[{f}, 
  f = zpq[p,q]/dpq[p, q] + 
    Sqrt[(zpq[p,q]/dpq[p, q]) + 1] Sqrt[(zpq[p,q]/dpq[p, q])  - 1]]
    
vpq[p_, q_] := 
 Module[{f}, 
  f = zpq[p, q]/dpq[p, q] + 
    Sqrt[(zpq[p, q]/dpq[p, q]) + 1] Sqrt[(zpq[p, q]/dpq[p, q])  - 1]]
    
    (*-----------------------------EXPANSION COEFFICIENTS---------------------------*)
    
    (*EtaFarPQ[n_, m_, p_, q_]*)
    
EtaFarPQ[n_, m_, p_, q_] := 
 Module[{f1, f2}, f2 = dd[p]/dd[q]; 
  f1 = If[p == q, 0, 
    n*Power[dd[p], n]*
     Power[-1, m] Sum[
      Power[dd[q] , 2 l + Abs[m]]* 
       Sum[((f2)^(2 k))/(Factorial[k]* Factorial[l - k]*Factorial[k + n]*Factorial[Abs[m] + l - k]), {k, 
         0, l}] (Factorial[n + Abs[m] + 2 l - 1]/
         Power[2 *zpq[p, q], n + Abs[m] + 2*l]), {l, 0, 90}]]]
         
(*MuFarPQ[n_, m_, p_, q_]*)
MuFarPQ[n_, m_, p_, q_] := 
 Module[{f1, f2, mm}, mm = Abs[m]; f2 = dd[p]/dd[q];
  f1 = If[p == q,0,
  	-2*n* 
     Power[dd[p], n] Power[-1, mm]  Sum[
      Power[dd[q], 2 l + mm]*
       Sum[Power[f2,2*k]/(Factorial[k]*Factorial[l - k]* 
           Factorial[k + n]*Factorial[Abs[m] + l - k]), {k, 0,l}] (Factorial[n + mm + 2 l]/Power[2* zpq[p, q], n + mm + 2*l + 1]), {l, 0, 90}]]]
           
(*EtaClosePQ[n_, m_, p_, q_]*)    
EtaClosePQ[n_, m_, p_, q_] := 
 Module[{f1, f2, d1, d2}, f2 = dd[p]/dd[q]; d1 = dd[p]/dpq[p, q]; 
  d2 = dd[q]/dpq[p, q]; 
  f1 = If[p == q, 0, 
    n *Power[d1, n]*Power[-1, m]*
     Sum[Power[vpq[p, q], -(n + m + 2*j)] *
       Sum[(Power[-1, j - l]/Gamma[j - l + 1])*(Power[d2, 
           2* l + m]) Sum[
          Power[f2, 
            2*k]/(Factorial[k]*Factorial[l - k]*Factorial[k + n]*
             Factorial[Abs[m] + l - k]), {k, 0, l}]   *
         Gamma[n + m + l + j ], {l, 0, j}], {j, 0, 19}] ]]
      
      
(*MuClosePQ[n_, m_, p_, q_]*)   
MuClosePQ[n_, m_, p_, q_] := 
 Module[{f1, f2, vpqprime}, f2 = dd[p]/dd[q]; 
  vpqprime = (1/(dpq[p, 
        q])) (1 +  ((zpq[p, q]/
          dpq[p, q])/(Sqrt[(zpq[p, q]/dpq[p, q]) + 
            1] Sqrt[(zpq[p, q]/dpq[p, q])  - 1])));
  f1 = If[p == q, 
    0, -n*vpqprime Power[dd[p]/dpq[p, q], n]*Power[-1, m]*
     Sum[ (n + m + 2*j)*Power[ vpq[p, q], -(n + m + 2*j + 1)]* 
       Sum[(Power[-1, j - l]/Gamma[j - l + 1]) Power[dd[q]/dpq[p, q], 
          2* l + m] Sum[
          Power[f2, 
            2*k]/(Factorial[k]*Factorial[l - k]*Factorial[k + n]*
             Factorial[Abs[m] + l - k]), {k, 0, l}] *  
              Gamma[n + m + l + j ], {l, 0, j}], {j, 0, 19}]]]   
         
vn[z_, n_, p_] := 
  Module[{f}, 
    f = Exp[-n ArcCosh[z/dd[p]]]]

EtaIntgPQ[n_, m_, p_, q_] := 
  Module[{f}, 
    f = If[p == q, 0, (1/Pi)*
      	NIntegrate[vn[dd[q] Cos[Etatemp] + zpq[p, q], n, p] Cos[m*Etatemp], {Etatemp,0, Pi},MaxRecursion->6]]]
        
MuIntgPQ[n_, m_, p_, q_] := 
  Module[{f, f1}, 
    f1 = vn[dd[q] Cos[Etatemp] + zpq[p, q], 1, p];
    f = If[p == q, 0, (2*n/dd[p])*(1/Pi)*
 	    NIntegrate[vn[dd[q] Cos[Etatemp] + zpq[p, q], n,p] Cos[m*Etatemp]/(vn[dd[q] Cos[Etatemp] + zpq[p, q], 1, p] - 
         1/vn[dd[q] Cos[Etatemp] + zpq[p, q], 1, p]), {Etatemp,0, Pi},MaxRecursion->6]]]
  



End[] 
Protect @@ Names["ExpansionCoefficient`*"]
EndPackage[] 

BeginPackage["ExpansionCoefficient`",{"geometry`"}] 
EndPackage[] 

BeginPackage["ExpansionCoefficient`",{"geometry`"}] 
EndPackage[] 








