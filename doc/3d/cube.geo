//+
SetFactory("OpenCASCADE");
Box(1) = {-10, -4, -4, 20, 8, 8};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Boundary") = {1, 4, 6, 5, 3, 2};
