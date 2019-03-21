//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};

//+
Physical Surface("domain") = {1};
//+
Physical Curve("topSurface") = {4, 3, 2, 1};
