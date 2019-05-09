Merge "test.msh";
Point(1) = {-5, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Curve("Domain") = {1};
//+
Physical Point("Boundary") = {1, 2};
