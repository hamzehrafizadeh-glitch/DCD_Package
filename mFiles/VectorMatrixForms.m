(* ::Package:: *)

BeginPackage["VectorMatrixForms`"]
Unprotect @@ Names["VectorMatrixForms`*"];
ClearAll @@ Names["VectorMatrixForms`*"];


VecCreate::usage=" VecCreate[funcname_,func_,n_,ntot_] VecCreate function is a general function that produces a vector by having func[n,q] and return vector named funcname which is a (n*ntot)x1 dimension."; 
MatrixCreate::usage=" MatrixCreate(matname,Rmnpq,n,ntot) This function generates a matrix of (n ntot)(n ntot) dimensions with components of R(m,n,p,q)."; 
VecRev::usage=" VecRev(funcname,func,n,ntot) This is the reverse of the vector function.It allocates vector terms func(n,q) according to the specific array of n*ntot length."; 



Begin["`Private`"] 


(*Change the path into mFiles path.*)


Get["geometry.m"];
Needs["geometry`"];



(*VecCreate function is a general function that produces a vector by having func[n,q] and return vector named "funcname" which is a (n*ntot)x1 dimension vector in the below form:
 (\[NoBreak]ReFn1[1,1]
ReFn1[2,1]
ReFn1[3,1]
ReFn1[1,2]
ReFn1[2,2]
ReFn1[3,2]

\[NoBreak]). Here we have 2 inclusions and n-the umber of harmonic terms; n is 3 here.*)
 
 VecCreate[funcname_,func_,n_,ntot_]:= Module[{f,ndim,i,j,ind1},
ndim = n*ntot;
ind1=0;
funcname=ConstantArray[0,{ndim}];
f=
For[i=1,i<ntot+1,i++,
     For[j=1,j<n+1,j ++,
ind1=j+n(i-1);
(*Print["i: ", i," j: ",j, " ind1: ",ind1];*)
funcname[[ind1]]=func[j,i]];
]]

(*Test it
Clear[funcname]
VecCreate[funcname,ReFn1,n,ntot]
MatrixForm[funcname]*)





(*This is the reverse of the vector function.It allocates vector terms (func[n,q]) according to the specific array of n*ntot length.*)
VecRev[funcname_,func_,n_,ntot_]:= Module[{f,ndim,i,j,ind1},
ndim = n*ntot;
ind1=0;
f=
For[i=1,i<ntot+1,i++,
     For[j=1,j<n+1,j ++,
ind1=j+n(i-1);
(*Print["i: ", i," j: ",j, " ind1: ",ind1];*)
func[j,i]=funcname[[ind1]]];
]]


(*MatrixCreate[matname_,Rmnpq_,n_,ntot_]
This function generates a matrix of (n*ntot)x(n*ntot) dimensions with components of R[m,n,p,q] in the below form.Here ntot=2 and n=2.
(\[NoBreak]Rmnpq[1,1,1,1]	Rmnpq[2,1,1,1]	Rmnpq[3,1,1,1]	Rmnpq[1,1,2,1]	Rmnpq[2,1,2,1]	Rmnpq[3,1,2,1]
Rmnpq[1,2,1,1]	Rmnpq[2,2,1,1]	Rmnpq[3,2,1,1]	Rmnpq[1,2,2,1]	Rmnpq[2,2,2,1]	Rmnpq[3,2,2,1]
Rmnpq[1,3,1,1]	Rmnpq[2,3,1,1]	Rmnpq[3,3,1,1]	Rmnpq[1,3,2,1]	Rmnpq[2,3,2,1]	Rmnpq[3,3,2,1]
Rmnpq[1,1,1,2]	Rmnpq[2,1,1,2]	Rmnpq[3,1,1,2]	Rmnpq[1,1,2,2]	Rmnpq[2,1,2,2]	Rmnpq[3,1,2,2]
Rmnpq[1,2,1,2]	Rmnpq[2,2,1,2]	Rmnpq[3,2,1,2]	Rmnpq[1,2,2,2]	Rmnpq[2,2,2,2]	Rmnpq[3,2,2,2]
Rmnpq[1,3,1,2]	Rmnpq[2,3,1,2]	Rmnpq[3,3,1,2]	Rmnpq[1,3,2,2]	Rmnpq[2,3,2,2]	Rmnpq[3,3,2,2]

\[NoBreak]).
Also,line If[p==q,Break[]] Clear[jnd]; guarantees if both terms are the same,the matrix element becomes zero.*)

Clear[MatrixCreate,matname]
MatrixCreate[matname_,Amnpq_,n_,ntot_]:=Module[{f, ndim,ndim2,i,j,p,q,m,ind,ind2,jnd},
ndim = 1*n*ntot;
     ndim2=ndim/2;
ind=2;
     m=0;
matname=ConstantArray[0,{ndim,ndim}];
For[q=1,q<ntot+1,q++,
     For[j=1,j<n+1,j ++,

For[p=1,p<=ntot,p++,For[i=1,i<n+1,i++,
(*If[p==q,Break[]] Clear[jnd];*)
jnd=n*(p-1)+(i-1)+1;
m+=1;
(*, ind-1,jnd*)
	matname[[ind-1,jnd]]=Amnpq[i,j,p,q]
];];
Clear[ind2];ind2=ind;
Clear[ind];ind=ind2+1;

];];]

(* Test it by
Clear[matname]
MatrixCreate[matname,Amnpq,2,3]
MatrixForm[matname]
Dimensions[matname]*)



End[] 
Protect @@ Names["VectorMatrixForms`*"]
EndPackage[] 

BeginPackage["VectorMatrixForms`",{"geometry`"}] 
EndPackage[] 

BeginPackage["VectorMatrixForms`",{"geometry`"}] 
EndPackage[] 











