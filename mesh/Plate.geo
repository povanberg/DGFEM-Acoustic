// Parameters

w = 5;        // Width of the plate
h = 5;        // Height of the playe
lc = 0.5;       // Catacteristic lenght
n = 0;        // Number of holes
r = 0.5;      // Radius of the holes

// Creation of the rectangle

Point(1) = {0,0,0,lc};
Point(2) = {w,0,0,lc};
Point(3) = {w,h,0,lc};
Point(4) = {0,h,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Macro to generate a hole

Macro Hcreate

  p1 = newp;  Point(newp) = {x,y,0,lc};
  p2 = newp;  Point(newp) = {x+r,y,0,lc};
  p3 = newp;  Point(newp) = {x,y+r,0,lc};
  p4 = newp;  Point(newp) = {x-r,y,0,lc};
  p5 = newp;  Point(newp) = {x,y-r,0,lc};

  c1 = newl;  Circle(c1) = {p2,p1,p3};
  c2 = newl;  Circle(c2) = {p3,p1,p4};
  c3 = newl;  Circle(c3) = {p4,p1,p5};
  c4 = newl;  Circle(c4) = {p5,p1,p2};

  Hcurves[] += {c1,c2,c3,c4};                       // Stores the curves of the hole in a list
  cl1 = newll; Curve Loop(cl1) = {c1,c2,c3,c4};     // Creates a curve loop of the hole
  Hloops[] += cl1;                                  // Stores the curve loop in a list

Return

// Creation of n holes

For t In {1:n:1}

  If (t<4)

    x = r+2*t*r;
    y = r+2*t*r;
    Call Hcreate;

  ElseIf (t==4)

    x = r+2*r;
    y = h-(r+2*r);
    Call Hcreate;

  Else

    x = r+6*r;
    y = h-(r+6*r);
    Call Hcreate;

  EndIf

EndFor

// Creation of the plate surface

cl1 = newll; Curve Loop(cl1) = {1,2,3,4};                           // Exterior curve loop of the plate
s1 = news;  Plane Surface(s1) = {cl1,Hloops[]};                     // Plane surface of the plate and removes the holes
ps1 = newreg; Physical Surface("Plate surface",ps1) = {s1};         // Plate physical surface

pc1 = newreg; Physical Curve("Holes border",pc1) = {Hcurves[]};     // Holes border physical curve
pc2 = newreg; Physical Curve("Infinity border",pc2) = {2,3,4};      // Infinity physical curve curve
pc2 = newreg; Physical Curve("Ground border",pc2) = {1};            // Ground physical curve

// Creates 3 physical curve groups for the holes

Hphase = Round(n/3);   // Size of a group

For t In {0:n-1:1}

  If (t<Hphase)
    Plist1[] += {Hcurves[t*4],Hcurves[t*4+1],Hcurves[t*4+2],Hcurves[t*4+3]};
  
  ElseIf (t>=Hphase && t<2*Hphase)
    Plist2[] += {Hcurves[t*4],Hcurves[t*4+1],Hcurves[t*4+2],Hcurves[t*4+3]};

  Else
    Plist3[] += {Hcurves[t*4],Hcurves[t*4+1],Hcurves[t*4+2],Hcurves[t*4+3]};
  EndIf

EndFor

Physical Curve("Phase 1",newreg) = {Plist1[]};      // First holes physical curve group
Physical Curve("Phase 2",newreg) = {Plist2[]};      // Second holes physical curve group
Physical Curve("Phase 3",newreg) = {Plist3[]};      // Third holes physical curve group

// Meshes the surface

Mesh ps1;