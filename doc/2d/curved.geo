//+
SetFactory("OpenCASCADE");
Circle(1) = {0, -0, 0, 5, 0, 2*Pi};
//+
Circle(2) = {2.5, -0, -0.1, 1.5, 0, 2*Pi};
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Surface(1) = {1, 2};
//+
Curve Loop(3) = {1};
//+
Physical Curve("Boundary") = {2, 1};
//+
Curve Loop(4) = {2};
//+
Curve Loop(5) = {1};
//+
Curve Loop(6) = {2};
//+
Plane Surface(1) = {5, 6};
//+
Physical Surface("Domain") = {1};

Mesh.SecondOrderLinear = 0;
