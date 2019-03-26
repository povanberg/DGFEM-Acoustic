Merge "circle.gmsh";
SetFactory("OpenCASCADE");
Circle(1) = {1.2, -1.3, 0, 5, 0, 2*Pi};
//+
Circle(2) = {3.7, -1.3, 0, 1.3, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Boundary") = {2, 1};
//+
Physical Surface("Domain") = {1};
