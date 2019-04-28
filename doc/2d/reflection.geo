//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Domain") = {1};
//+
Physical Curve("Reflecting") = {1};
