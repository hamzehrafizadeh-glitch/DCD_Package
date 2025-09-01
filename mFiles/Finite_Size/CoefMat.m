(* ::Package:: *)

 BeginPackage[ "CoefMat`"];
 Unprotect@@Names["CoefMat`*"];
 ClearAll@@Names["CoefMat`*"];


(*Uij::usage ="Sxx Stress Field.";
Cij::usage ="Syy Stress Field.";*)
dxdy::usage ="It determines dr=(dx,dy) between two elements i and j.";
dx::usage ="It determines dx term of dr=(dx,dy).";
dy::usage ="It determines dy term of dr=(dx,dy).";



 Begin[ "Private`"];
 


(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["geometry.m"];
Needs["geometry`"];
Get["FieldsBE.m"];
Needs["FieldsBE`"];
XY =<< XY.dat;


(*Change m to Mu*)
(*Chi=3-4*Nu;Nu=0.07;E1=100.1;(*m=E1/(2*(Nu+1))*)
Mu=570;*)


k=Chi0;m=Mu0;   

(*dx between two elements on the boundary.*)
dxdy[i_,j_]:=Block[{f},f=Rotate2Dcc[(XY[[i,1]]-XY[[j,1]]),(XY[[i,2]]-XY[[j,2]]),XY[[j,4]]]]

dx[i_,j_]:=Re[dxdy[i,j]]
dy[i_,j_]:=Im[dxdy[i,j]]

       
 
Print["Matrixes Uij and Cij Start."];
Print["Time"];
Print[DateList[]];
Uij=Table[0,{i,2*NTOTAL},{j,1,2*NTOTAL}];
Cij=Table[0,{i,2*NTOTAL},{j,1,2*NTOTAL}];

For[i=1,i<NTOTAL+1,i++,
For[j=1,j<NTOTAL+1,j++,      
cr=If [dy[i,j]==0,10^-13,0];

cr=If [dy[i,j]==0,-10^-13,0];
 sx1=SxBE[dx[i,j],dy[i,j]+cr,1,0,XY[[j,3]],k];
        sx2=SxBE[dx[i,j],dy[i,j]+cr,0,1,XY[[j,3]],k];
        sy1=SyBE[dx[i,j],dy[i,j]+cr,1,0,XY[[j,3]],k];
        sy2=SyBE[dx[i,j],dy[i,j]+cr,0,1,XY[[j,3]],k];
        txy1=TxyBE[dx[i,j],dy[i,j]+cr,1,0,XY[[j,3]],k];
        txy2=TxyBE[dx[i,j],dy[i,j]+cr,0,1,XY[[j,3]],k];
        
        ux1=UxBE[dx[i,j],dy[i,j]+cr,1,0,XY[[j,3]]*0.99999,k,m];
        ux2=UxBE[dx[i,j],dy[i,j]+cr,0,1,XY[[j,3]]*0.99999,k,m];
        uy1=UyBE[dx[i,j],dy[i,j]+cr,1,0,XY[[j,3]]*0.99999,k,m];
        uy2=UyBE[dx[i,j],dy[i,j]+cr,0,1,XY[[j,3]]*0.99999,k,m];
        
(* sqw=-XY[[j,4]]+XY[[i,4]];*)
        
        bthi=-XY[[j,4]];
        athi=XY[[i,4]] ;  sqw=athi+bthi;

 Cij[[i,j]]=(sx1*Sin[sqw]^2+sy1*Cos[sqw]^2-2*txy1*Sin[sqw]*Cos[sqw]);
        (*s22_1 rotate*)
        Cij[[i+NTOTAL,j]]=-(sx1-sy1)*Sin[sqw]*Cos[sqw]+txy1*(Cos[sqw]^2-Sin[sqw]^2);
     (*s12_1 rotate*)
        Cij[[i,j+NTOTAL]]=(sx2*Sin[sqw]^2+sy2*Cos[sqw]^2-2*txy2*Sin[sqw]*Cos[sqw]);
         (*s22_2 rotate*)
        Cij[[i+NTOTAL,j+NTOTAL]]=-(sx2-sy2)*Sin[sqw]*Cos[sqw]+txy2*(Cos[sqw]^2-Sin[sqw]^2);
   (*s12_2 rotate*)
        

(*It is also required to rotate displacement tensor,but in its local coordinate system.*)
Uij[[i,j]]=Re[Rotate2Dcc[ux1,uy1,-XY[[j,4]]]];
Uij[[i+NTOTAL,j]]=Im[Rotate2Dcc[ux1,uy1,-XY[[j,4]]]];

Uij[[i,j+NTOTAL]]=Re[Rotate2Dcc[ux2,uy2,-XY[[j,4]]]];
Uij[[i+NTOTAL,j+NTOTAL]]=Im[Rotate2Dcc[ux2,uy2,-XY[[j,4]]]];
]]
       


Do[If[Abs[Cij[[i,j]]]<=10^-12,Cij[[i,j]]=0],{i,1,2*NTOTAL},{j,1,2*NTOTAL}]
Do[If[Abs[Uij[[i,j]]]<=10^-12,Uij[[i,j]]=0],{i,1,2*NTOTAL},{j,1,2*NTOTAL}]



Uij >> Uij.dat;
Cij >> Cij.dat;

Print["Matrixes Uij and Cij are built."];
Print["Time"];
Print[DateList[]];


 End[];
 Protect@@Names["CoefMat`*"]
 EndPackage[]

BeginPackage["CoefMat`",{"GenerateMesh`","FieldsBE`","geometry`"}] 
EndPackage[] 








































