Include "2DMacro.geo";
Include "2DParam.pro";

/* Borders */

pBord[0] = newp; Point(pBord[0]) = {-xBord/2,-yBord/2,0,LcBord};
pBord[1] = newp; Point(pBord[1]) = {xBord/2,-yBord/2,0,LcBord};
pBord[2] = newp; Point(pBord[2]) = {xBord/2,yBord/2,0,LcBord};
pBord[3] = newp; Point(pBord[3]) = {-xBord/2,yBord/2,0,LcBord};
pBord[4] = pBord[0];

For i In{0:3}
    lBord[i] = newl; Line(lBord[i]) = {pBord[i],pBord[i+1]};
EndFor

/* Holes */

clHoles = {};
cHoles = {};

For i In {1:nHoleCircle}
    x = XC~{i};
    y = YC~{i};
    z = ZC~{i};
    r = RC~{i};
    Call "HoleCircle";
EndFor

For i In {1:nHoleLong}
    x = XL~{i};
    y = YL~{i};
    z = ZL~{i};
    r = RL~{i};
    w = WL~{i};
    Call "HoleLong";
EndFor

clBord = newll; Curve Loop(clBord) = {lBord[]};
sBord = news; Plane Surface(sBord) = {clBord,clHoles[]};

/* Physical Entities */

Physical Surface("Domain",1000) = {sBord};

If(BcBord==1 && BcHole==1)
    Physical Curve("Absorbing",2000) = {lBord[],cHoles[]};
ElseIf(BcBord==0 && BcHole==0)
    Physical Curve("Reflecting",2000) = {lBord[],cHoles[]};
ElseIf(BcBord==1 && BcHole==0)
    Physical Curve("Absorbing",2000) = {lBord[]};
    Physical Curve("Reflecting",3000) = {cHoles[]};
ElseIf(BcBord==0 && BcHole==1)
    Physical Curve("Absorbing",2000) = {cHoles[]};
    Physical Curve("Reflecting",3000) = {lBord[]};
EndIf
