//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-2.5, -2.5, 0, 5, 5, 0};
//+
Circle(5) = {1.2, 1.2, 0, 0.5, 0, 2*Pi};

//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Curve Loop(3) = {5};
//+
Plane Surface(2) = {2, 3};
//+
Physical Surface("Domain") = {2};
//+
Physical Curve("Reflecting") = {4, 5, 1, 2, 3};
