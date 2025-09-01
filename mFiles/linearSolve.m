(* ::Package:: *)

BeginPackage["linearSolve`"]
Unprotect@@Names["linearSolve`*"];
ClearAll@@Names["linearSolve`*"];



abmnpq`a::usage=" a"; 
abmnpq`b::usage=" b"; 
abmnpq`c::usage=" c"; 
abmnpq`d::usage=" d"; 
xbartemp::usage=" d"; 

Begin["`Private`"] 

(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
(*Get["TemperatureInclusions.m"];
Needs["TemperatureInclusions`"];*)
(*Get["farfield.m"];*)
Needs["matrix`"]
Needs["farfield`"];


Unprotect[a,b,c,d]
(*Print["abar: ", abar]*)
PutAppend[ LinearAlgebra`MatrixConditionNumber[abar],FileNameJoin[{FilePath, "output", "matrix"}]]
xbartemp=LinearSolve[abar,bcbar]

For[j=1,j<ntot+1,j++,For[i=1,i<n+1,i++,a1[i,j]=xbartemp[[1+4 (i-1)+(j-1) (4 n)]];
b1[i,j]=xbartemp[[2+4 (i-1)+(j-1) (4 n)]];
c1[i,j]=xbartemp[[3+4 (i-1)+(j-1) (4 n)]];
d1[i,j]=xbartemp[[4+4 (i-1)+(j-1) (4 n)]];
a2[i,j]=xbartemp[[1+4 (i-1)+(j-1+ntot) (4 n)]];
b2[i,j]=xbartemp[[2+4 (i-1)+(j-1+ntot) (4 n)]];
c2[i,j]=xbartemp[[3+4 (i-1)+(j-1+ntot) (4 n)]];
d2[i,j]=xbartemp[[4+4 (i-1)+(j-1+ntot) (4 n)]];
abmnpq`a[i,j]=a1[i,j]+I a2[i,j];abmnpq`b[i,j]=b1[i,j]+I b2[i,j];
abmnpq`c[i,j]=c1[i,j]+I c2[i,j];abmnpq`d[i,j]=d1[i,j]+I d2[i,j]]];

Print["LinearSolve.m"];


Clear[i,j]
s=OpenWrite[FileNameJoin[{FilePath, "output", "linearSolve"}]]
For[j=1,j<ntot+1,j++,For[i=1,i<n+1,i++,
Write[s,abmnpq`a[i,j]];
Write[s,abmnpq`b[i,j]];
Write[s,abmnpq`c[i,j]];
Write[s,abmnpq`d[i,j]];
]]
Close[s]



(*Print[i]
For[j = 1, j < 5 + 1, j++, 
 For[i = 1, i < 5 + 1, i++, 
  PutAppend[{o1[i 5, j], o2[i , j], o3[i , j], o4[i , j]}, 
   "~/Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/\
Package/linearSolve/test"] ]]*)


End[] 
Protect@@Names["linearSolve`*"]
EndPackage[] 

BeginPackage["linearSolve`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`" ,"VectorMatrixForms`", "matrix`","TemperatureInclusions`" , "farfield`"}] 
EndPackage[] 
























