//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Surface("Domain") = {1};
//+
Physical Curve("topSurface") = {3};
//+
Physical Curve("rightSurface") = {2};
//+
Physical Curve("bottomSurface") = {1};
//+
Physical Curve("leftSurface") = {4};
