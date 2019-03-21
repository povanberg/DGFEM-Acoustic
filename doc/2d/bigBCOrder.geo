//+
SetFactory("OpenCASCADE");
Rectangle(1) = {5, -5, 0, 10, 10, 0};
//+
Physical Curve("leftSurface") = {4};
//+
Physical Curve("topSurface") = {3};
//+
Physical Curve("bottomSurface") = {1};
//+
Physical Curve("rightSurface") = {2};
//+
Physical Surface("Domain") = {1};
