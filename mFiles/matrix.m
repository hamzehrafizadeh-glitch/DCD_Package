(* ::Package:: *)

BeginPackage["matrix`"]
Unprotect @@ Names["matrix`*"];
ClearAll @@ Names["matrix`*"];

abar::usage=" The problem is reduced to a set of linear system of equations. To solve the system using linear solver, it is written in matrix format. "; 

Begin["`Private`"] 

(*Print["Path", $Path];*)
Get["abmnpq.m"]
Needs["abmnpq`"];

alpha11[m_,n_,p_,q_]:=Block[{f,f1},f1=v0[q]^2-v0[q]^(-2);f=-Chi0*smnpq[m,n,p,q]+Power[v0[q],2*n]*(Conjugate[bmnpqa[m,n,p,q]-n*f1*smnpq[m,n,p,q]])]
alpha12[m_,n_,p_,q_]:=Block[{f},f=Power[v0[q],2*n]*Conjugate[bmnpqb[m,n,p,q]]]
alpha13[m_,n_,p_,q_]:=Block[{f,f1},f1=v0[q]^2-v0[q]^(-2);f=I*(-Chi0*smnpq[m,n,p,q]-Power[v0[q],2*n]*(Conjugate[bmnpqa[m,n,p,q]-n*f1*smnpq[m,n,p,q]]))]
alpha14[m_,n_,p_,q_]:=Block[{f},f=-I Power[v0[q],2*n]*Conjugate[bmnpqb[m,n,p,q]]]


alpha21[m_,n_,p_,q_]:=Block[{f},f=Chi0*Power[v0[q],2*n]*smnpq[m,n,p,q]-Conjugate[bmnpqa[m,n,p,q]]]
alpha22[m_,n_,p_,q_]:=Block[{f},f=-Conjugate[bmnpqb[m,n,p,q]]]
alpha23[m_,n_,p_,q_]:=Block[{f},f=I*(f=Chi0*Power[v0[q],2*n]*smnpq[m,n,p,q]+Conjugate[bmnpqa[m,n,p,q]])]
alpha24[m_,n_,p_,q_]:=Block[{f},f=I Conjugate[bmnpqb[m,n,p,q]]]


alpha31[m_,n_,p_,q_]:=Block[{f,f1},f1=v0[q]^2-v0[q]^(-2);f=-smnpq[m,n,p,q]-Power[v0[q],2*n]*(Conjugate[bmnpqa[m,n,p,q]-n*f1*smnpq[m,n,p,q]])]
alpha32[m_,n_,p_,q_]:=Block[{f},f=-Power[v0[q],2*n]*Conjugate[bmnpqb[m,n,p,q]]]
alpha33[m_,n_,p_,q_]:=Block[{f,f1},f1=v0[q]^2-v0[q]^(-2);f=I*(-smnpq[m,n,p,q]+Power[v0[q],2*n]*(Conjugate[bmnpqa[m,n,p,q]-n*f1*smnpq[m,n,p,q]]))]
alpha34[m_,n_,p_,q_]:=Block[{f},f=I*Power[v0[q],2*n]*Conjugate[bmnpqb[m,n,p,q]]]


alpha41[m_,n_,p_,q_]:=Block[{f},f=-Power[v0[q],2*n]*smnpq[m,n,p,q]-Conjugate[bmnpqa[m,n,p,q]]]
alpha42[m_,n_,p_,q_]:=Block[{f},f= -Conjugate[bmnpqb[m,n,p,q]]]
alpha43[m_,n_,p_,q_]:=Block[{f},f=I(-Power[v0[q],2*n]*smnpq[m,n,p,q]+Conjugate[bmnpqa[m,n,p,q]])]
alpha44[m_,n_,p_,q_]:=Block[{f},f= I Conjugate[bmnpqb[m,n,p,q]]]




Clear[ll0]
ndim = 8*n*ntot;
ndim2=ndim/2;
ll0=ConstantArray[0,{ndim,ndim}];
Clear[i,j]
ind=1;
For[q=1,q<ntot+1,q++,
     For[j=1,j<n+1,j ++,

(*"Matrix for real part"*)
(*Print[j];*)
ll0[[ind,ind]]=Re[Chi0];
ll0[[ind,ind+1]]=0;
ll0[[ind,ind+2]]=Re[-2*j* Sinh[2 Zeta0[q]]*Power[v0[q],2*j]-Chi[q]];
ll0[[ind,ind+3]]=Re[Power[v0[q],2*j]];

ll0[[ind+1,ind]]=0;
ll0[[ind+1,ind+1]]=1;
ll0[[ind+1,ind+2]]=Re[Power[v0[q],2*j]* Chi[q]];
ll0[[ind+1,ind+3]]=-1;

ll0[[ind+2,ind]]=1;
ll0[[ind+2,ind+1]]=0;
ll0[[ind+2,ind+2]]=Re[-MuBar[q]+2*j* Sinh[2 Zeta0[q]] Power[v0[q],2*j]* MuBar[q]];
ll0[[ind+2,ind+3]]=-Re[Power[v0[q],2*j]* MuBar[q]];

ll0[[ind+3,ind]]=0;
ll0[[ind+3,ind+1]]=1;
ll0[[ind+3,ind+2]]=-Re[Power[v0[q],2*j]* MuBar[q]];
ll0[[ind+3,ind+3]]=-Re[MuBar[q]];


ind+=ntot*n*4;


(*"Matrix for right hand "*)

  ll0[[ind,ind]]=Re[Chi0];
ll0[[ind,ind+1]]=0;
ll0[[ind,ind+2]]=Re[2*j* Sinh[2 Zeta0[q]] Power[v0[q],2*j]-Chi[q]];
ll0[[ind,ind+3]]=-Re[Power[v0[q],2*j]];

ll0[[ind+1,ind]]=0;
ll0[[ind+1,ind+1]]=-1;
ll0[[ind+1,ind+2]]=Re[Power[v0[q],2*j]* Chi[q]];
ll0[[ind+1,ind+3]]=1;

ll0[[ind+2,ind]]=1;
ll0[[ind+2,ind+1]]=0;
ll0[[ind+2,ind+2]]=Re[-MuBar[q]-2 *j*Sinh[2 Zeta0[q]]* Power[v0[q],2*j]* MuBar[q]];
ll0[[ind+2,ind+3]]=Re[Power[v0[q],2*j]* MuBar[q]];

ll0[[ind+3,ind]]=0;
ll0[[ind+3,ind+1]]=-1;
ll0[[ind+3,ind+2]]=-Re[Power[v0[q],2*j]* MuBar[q]];
ll0[[ind+3,ind+3]]=Re[MuBar[q]];

ind-=ntot*n*4-4;

]]

Clear[ll]
ll=ConstantArray[0,{ndim,ndim}];
Clear[i,j]
ind=5;
For[q=1,q<ntot+1,q++,
     For[j=1,j<n+1,j ++,

For[p=1,p<=ntot,p++,For[i=1,i<n+1,i++,If[p==q,Break[]] Clear[jnd];jnd=4*n*(p-1)+4*(i-1)+1;

ll[[ind-4,jnd]]=Re[alpha11[i,j,p,q]];
ll[[ind-4,jnd+1]]=Re[alpha12[i,j,p,q]];
ll[[ind-4+ndim2,jnd]]=Im[alpha11[i,j,p,q]];
ll[[ind-4+ndim2,jnd+1]]=Im[alpha12[i,j,p,q]];
ll[[ind-3,jnd]]=Re[alpha21[i,j,p,q]];
ll[[ind-3,jnd+1]]=Re[alpha22[i,j,p,q]];
ll[[ind-3+ndim2,jnd]]=Im[alpha21[i,j,p,q]];
ll[[ind-3+ndim2,jnd+1]]=Im[alpha22[i,j,p,q]];
ll[[ind-2,jnd]]=Re[alpha31[i,j,p,q]];
ll[[ind-2,jnd+1]]=Re[alpha32[i,j,p,q]];
ll[[ind-2+ndim2,jnd]]=Im[alpha31[i,j,p,q]];
ll[[ind-2+ndim2,jnd+1]]=Im[alpha32[i,j,p,q]];
ll[[ind-1,jnd]]=Re[alpha41[i,j,p,q]];
ll[[ind-1,jnd+1]]=Re[alpha42[i,j,p,q]];
ll[[ind-1+ndim2,jnd]]=Im[alpha41[i,j,p,q]];
ll[[ind-1+ndim2,jnd+1]]=Im[alpha42[i,j,p,q]];


ll[[ind-4,jnd+ndim2]]=Re[alpha13[i,j,p,q]];
ll[[ind-4,jnd+1+ndim2]]=Re[alpha14[i,j,p,q]];
ll[[ind-4+ndim2,jnd+ndim2]]=Im[alpha13[i,j,p,q]];
ll[[ind-4+ndim2,jnd+1+ndim2]]=Im[alpha14[i,j,p,q]];
ll[[ind-3,jnd+ndim2]]=Re[alpha23[i,j,p,q]];
ll[[ind-3,jnd+1+ndim2]]=Re[alpha24[i,j,p,q]];
ll[[ind-3+ndim2,jnd+ndim2]]=Im[alpha23[i,j,p,q]];
ll[[ind-3+ndim2,jnd+1+ndim2]]=Im[alpha24[i,j,p,q]];
ll[[ind-2,jnd+ndim2]]=Re[alpha33[i,j,p,q]];
ll[[ind-2,jnd+1+ndim2]]=Re[alpha34[i,j,p,q]];
ll[[ind-2+ndim2,jnd+ndim2]]=Im[alpha33[i,j,p,q]];
ll[[ind-2+ndim2,jnd+1+ndim2]]=Im[alpha34[i,j,p,q]];
ll[[ind-1,jnd+ndim2]]=Re[alpha43[i,j,p,q]];
ll[[ind-1,jnd+1+ndim2]]=Re[alpha44[i,j,p,q]];
ll[[ind-1+ndim2,jnd+ndim2]]=Im[alpha43[i,j,p,q]];
ll[[ind-1+ndim2,jnd+1+ndim2]]=Im[alpha44[i,j,p,q]];]]



Clear[ind2];ind2=ind;
Clear[ind];ind=ind2+4;

]]
abar=-ll+ll0;
End[] 




Protect@@Names["matrix`*"]
EndPackage[] 

BeginPackage["matrix`",{"geometry`", "ExpansionCoefficient`" , "abmnpq`"}] 
EndPackage[] 







