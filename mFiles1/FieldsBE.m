(* ::Package:: *)

 BeginPackage[ "FieldsBE`"];
 Unprotect@@Names["FieldsBE`*"];
 ClearAll@@Names["FieldsBE`*"];


SxBE::usage ="Sxx Stress Field.";
SyBE::usage ="Syy Stress Field.";
TxyBE::usage ="Txy Stress Field.";
UxBE::usage ="Ux Stress Field.";
UyBE::usage ="Uy Stress Field.";
Rotate2Dcc::usage ="2D rotation matrix. This rotates a vector from global coordinate system to the local one.";


 Begin[ "Private`"];
 
(*AppendTo[$Path,ToFileName[{$HomeDirectory,"Dropbox/Mathematica/ObliqueFiniteArray/Fracture/1Inclusion/Package/mFiles"}]];*)
Get["GenerateMesh.m"];
Needs["GenerateMesh`"];



(*m is a Shear modulus Mu;
k = \[Chi] or Chi; and a is a distance between any pair elements.
Po and Qo determine applied tenstion and shear stress respectively.*)


SxBE[x_,y_,Po_,Qo_,a_,k_]:=Block[{f},f=(1/(4*a*(1+k)*Pi))*(-2*(3+k)*Qo*(a-x)-4*(3+k)*Qo*x+2*(3+k)*Qo*(a+x)+8*Po*y-(4*y*((-a)*((-Po)*x+Qo*y)+Po*(x^2+y^2)))/(x^2+y^2)-(4*y*(a*((-Po)*x+Qo*y)+Po*(x^2+y^2)))/(x^2+y^2)+2*(a*(-3+k)*Po-(-3+k)*Po*x+(5+k)*Qo*y)*ArcTan[((a-x)/(0+y))]-2*(a*(-3+k)*Po+(-3+k)*Po*x-(5+k)*Qo*y)*ArcTan[(x/(0+y))]+2*(a*(-3+k)*Po-(-3+k)*Po*x+(5+k)*Qo*y)*ArcTan[(x/(0+y))]+2*(a*(-3+k)*Po+(-3+k)*Po*x-(5+k)*Qo*y)*ArcTan[((a+x)/(0+y))]+(a*(3+k)*Qo-(3+k)*Qo*x-(-5+k)*Po*y)*Log[((a-x)^2+y^2)]-(a*(3+k)*Qo-(3+k)*Qo*x-(-5+k)*Po*y)*Log[(x^2+y^2)]+(a*(3+k)*Qo+(3+k)*Qo*x+(-5+k)*Po*y)*Log[(x^2+y^2)]-(a*(3+k)*Qo+(3+k)*Qo*x+(-5+k)*Po*y)*Log[((a+x)^2+y^2)])];


SyBE[x_,y_,Po_,Qo_,a_,k_]:=Block[{f},f=(1/(4*a*(1+k)*Pi))*(2*(-1+k)*Qo*(a-x)+4*(-1+k)*Qo*x-2*(-1+k)*Qo*(a+x)-8*Po*y+(4*y*((-a)*((-Po)*x+Qo*y)+Po*(x^2+y^2)))/(x^2+y^2)+(4*y*(a*((-Po)*x+Qo*y)+Po*(x^2+y^2)))/(x^2+y^2)-2*(1+k)*(a*Po-Po*x+Qo*y)*ArcTan[((a-x)/y)]+2*(1+k)*(a*Po+Po*x-Qo*y)*ArcTan[(x/(0+y))]-2*(1+k)*(a*Po-Po*x+Qo*y)*ArcTan[(x/(0+y))]-2*(1+k)*(a*Po+Po*x-Qo*y)*ArcTan[((a+x)/(0+y))]-(-1+k)*(a*Qo-Qo*x-Po*y)*Log[((a-x)^2+y^2)]+(-1+k)*(a*Qo-Qo*x-Po*y)*Log[(x^2+y^2)]-(-1+k)*(a*Qo+Qo*x+Po*y)*Log[(x^2+y^2)]+(-1+k)*(a*Qo+Qo*x+Po*y)*Log[((a+x)^2+y^2)])];
TxyBE[x_,y_,Po_,Qo_,a_,k_]:=Block[{f},f=(-(1/(4*a*(1+k)*Pi)))*(2*(-1+k)*Po*(a-x)+4*(-1+k)*Po*x-2*(-1+k)*Po*(a+x)-8*Qo*y+(4*y*((-a)*(Qo*x+Po*y)+Qo*(x^2+y^2)))/(x^2+y^2)+(4*y*(a*(Qo*x+Po*y)+Qo*(x^2+y^2)))/(x^2+y^2)+2*(a*(1+k)*Qo-(1+k)*Qo*x-(-3+k)*Po*y)*ArcTan[((a-x)/(0+y))]+2*(a*(1+k)*Qo-(1+k)*Qo*x-(-3+k)*Po*y)*ArcTan[(x/(0+y))]-2*(a*(1+k)*Qo+(1+k)*Qo*x+(-3+k)*Po*y)*ArcTan[(x/(0+y))]+2*(a*(1+k)*Qo+(1+k)*Qo*x+(-3+k)*Po*y)*ArcTan[((a+x)/(0+y))]-(a*(-1+k)*Po-(-1+k)*Po*x+(3+k)*Qo*y)*Log[((a-x)^2+y^2)]-(a*(-1+k)*Po+(-1+k)*Po*x-(3+k)*Qo*y)*Log[(x^2+y^2)]+(a*(-1+k)*Po-(-1+k)*Po*x+(3+k)*Qo*y)*Log[(x^2+y^2)]+(a*(-1+k)*Po+(-1+k)*Po*x-(3+k)*Qo*y)*Log[((a+x)^2+y^2)])];





UxBE[x_,y_,Po_,Qo_,a_,k_,m_]:=Block[{f},f=1/(8 a (1+k) m \[Pi]) (-a^2 (1+k) Qo+2 a (a (Qo+2 k Qo)-k Qo x+2 Po y)-4 y (a (1+k) Qo-(1+k) Qo x+Po y) ArcTan[(a-x)/(0+y)]-4 y (a (1+k) Qo-(1+k) Qo x+Po y) ArcTan[x/(0+y)]-a^2 k Qo Log[(a-x)^2+y^2]-(2 a (k Qo x-Po y)+2 y (Po x+Qo y)+k Qo (-x^2+y^2)) Log[x^2+y^2]+(2 a (k Qo x-Po y)+2 y (Po x+Qo y)+k Qo (-x^2+y^2)) Log[a^2-2 a x+x^2+y^2])-1/(8 a (1+k) m \[Pi]) (a^2 (1+k) Qo-2 a (a (Qo+2 k Qo)+k Qo x-2 Po y)-4 y (a (1+k) Qo+(1+k) Qo x-Po y) ArcTan[x/(0+y)]+4 y (a (1+k) Qo+(1+k) Qo x-Po y) ArcTan[(a+x)/(0+y)]-(2 a (k Qo x-Po y)-2 y (Po x+Qo y)+k Qo (x^2-y^2)) Log[x^2+y^2]+(2 a (k Qo x-Po y)-2 y (Po x+Qo y)+k Qo (x^2-y^2)) Log[a^2+2 a x+x^2+y^2]+a^2 k Qo Log[(a+x)^2+y^2])];
UyBE[x_,y_,Po_,Qo_,a_,k_,m_]:=Block[{f},f=1/(8 a (1+k) m \[Pi]) (-a^2 (-1+k) Po-2 a (a (Po-2 k Po)+k Po x-2 Qo y)-4 y (a (-1+k) Po-(-1+k) Po x+Qo y) ArcTan[(a-x)/(0+y)]-4 y (a (-1+k) Po-(-1+k) Po x+Qo y) ArcTan[x/(0+y)]-a^2 k Po Log[(a-x)^2+y^2]+(2 y (-Qo x+Po y)+a (-2 k Po x+2 Qo y)+k Po (x^2-y^2)) Log[x^2+y^2]-(2 y (-Qo x+Po y)+a (-2 k Po x+2 Qo y)+k Po (x^2-y^2)) Log[a^2-2 a x+x^2+y^2])-1/(8 a (1+k) m \[Pi]) (a^2 (-1+k) Po-2 a (a (-1+2 k) Po+k Po x-2 Qo y)+4 y (a (Po-k Po)-(-1+k) Po x+Qo y) ArcTan[x/(0+y)]-4 y (a (Po-k Po)-(-1+k) Po x+Qo y) ArcTan[(a+x)/(0+y)]+(-2 y (-Qo x+Po y)+a (-2 k Po x+2 Qo y)+k Po (-x^2+y^2)) Log[x^2+y^2]-(-2 y (-Qo x+Po y)+a (-2 k Po x+2 Qo y)+k Po (-x^2+y^2)) Log[a^2+2 a x+x^2+y^2]+a^2 k Po Log[(a+x)^2+y^2])];


Rotate2Dcc[x_,y_,thi_]:=Block[{f,xp,yp},xp=x*Cos[thi]+y*Sin[thi];
yp=y*Cos[thi]-x*Sin[thi];f=xp+I*yp]



 End[];
 Protect@@Names["FieldsBE`*"]
 EndPackage[]

BeginPackage["FieldsBE`",{"GenerateMesh`"}] 
EndPackage[] 


































