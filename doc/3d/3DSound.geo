Bc = DefineNumber[0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "Parameters/1Boundary"];
Lc = DefineNumber[1.0, Name "Parameters/2Mesh size", Highlight "Green", Units "m", Min 1e-2, Max 10];
x = DefineNumber[10.0, Name "Parameters/3Width", Highlight "Green", Units "m", Min 0, Max 100];
y = DefineNumber[10.0, Name "Parameters/4Height", Highlight "Green", Units "m", Min 0, Max 100];
z = DefineNumber[10.0, Name "Parameters/4Thickness", Highlight "Green", Units "m", Min 0, Max 100];

pBoud[0] = newp; Point(pBoud[0]) = {-x/2,-y/2,-z/2,Lc};
pBoud[1] = newp; Point(pBoud[1]) = {x/2,-y/2,-z/2,Lc};
pBoud[2] = newp; Point(pBoud[2]) = {x/2,y/2,-z/2,Lc};
pBoud[3] = newp; Point(pBoud[3]) = {-x/2,y/2,-z/2,Lc};
pBoud[4] = pBoud[0];

pBoud[5] = newp; Point(pBoud[5]) = {-x/2,-y/2,z/2,Lc};
pBoud[6] = newp; Point(pBoud[6]) = {x/2,-y/2,z/2,Lc};
pBoud[7] = newp; Point(pBoud[7]) = {x/2,y/2,z/2,Lc};
pBoud[8] = newp; Point(pBoud[8]) = {-x/2,y/2,z/2,Lc};
pBoud[9] = pBoud[5];

For i In{0:3}
    lBoud[i] = newl; Line(lBoud[i]) = {pBoud[i],pBoud[i+1]};
    lBoud[i+5] = newl; Line(lBoud[i+5]) = {pBoud[i+5],pBoud[i+6]};
    lBoud[i+10] = newl; Line(lBoud[i+10]) = {pBoud[i],pBoud[i+5]};
EndFor

clBoud[0] = newll; Curve Loop(clBoud[0]) = {lBoud[0],lBoud[1],lBoud[2],lBoud[3]};
clBoud[1] = newll; Curve Loop(clBoud[1]) = {lBoud[5],lBoud[6],lBoud[7],lBoud[8]};
clBoud[2] = newll; Curve Loop(clBoud[2]) = {lBoud[10],lBoud[5],-lBoud[11],-lBoud[0]};
clBoud[3] = newll; Curve Loop(clBoud[3]) = {lBoud[11],lBoud[6],-lBoud[12],-lBoud[1]};
clBoud[4] = newll; Curve Loop(clBoud[4]) = {lBoud[12],lBoud[7],-lBoud[13],-lBoud[2]};
clBoud[5] = newll; Curve Loop(clBoud[5]) = {lBoud[13],lBoud[8],-lBoud[10],-lBoud[3]};

For i In{0:5}
    sBoud[i] = news; Plane Surface(sBoud[i]) = {clBoud[i]};
EndFor

slBoud = newsl; Surface Loop(slBoud) = {sBoud[]};
vBoud = newv; Volume(vBoud) = {slBoud};

If(Bc==0)
Physical Surface("Reflecting",1000) = {sBoud[]};
Else
    Physical Surface("Absorbing",1000) = {sBoud[]};
EndIf
Physical Volume("Domain",2000) = {vBoud};