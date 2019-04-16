//+
SetFactory("OpenCASCADE");
Rectangle(2) = {0, 0, 0, 20, 10, 0};
//+
Circle(5) = {10, 5, 0, 1, 0, 2*Pi};
//+
Circle(6) = {15, 7.5, 0, 1, 0, 2*Pi};
//+
Circle(7) = {15, 2.5, 0, 1, 0, 2*Pi};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Curve Loop(3) = {5};
//+
Curve Loop(4) = {6};
//+
Curve Loop(5) = {7};
//+
Plane Surface(1) = {2, 3, 4, 5};
//+
Physical Surface("Domain") = {1};
//+
Physical Curve("Absorbing") = {1, 4, 3, 2};
//+
Physical Curve("Reflecting") = {5, 6, 7};

