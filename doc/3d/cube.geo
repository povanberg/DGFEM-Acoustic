//+
SetFactory("OpenCASCADE");
Box(1) = {-10, -10, -10, 20, 20, 20};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Boundary") = {1, 4, 6, 5, 3, 2};
