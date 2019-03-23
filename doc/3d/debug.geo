//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Volume("domain") = {1};
//+
Physical Surface("Boundary") = {6, 2, 4, 5, 1, 3};
