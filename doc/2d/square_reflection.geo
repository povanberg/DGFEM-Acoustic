//+
SetFactory("OpenCASCADE");
out = news;
Rectangle(out) = {-5, -5, 0, 10, 10, 0};
//+
in = news;
Rectangle(in) = {-1, 2.5, 0, 2, 2, 0};
//+
domain = news;
BooleanDifference(domain) = { Surface{out}; Delete; }{ Surface{in}; Delete; };
//+
Physical Curve("Absorbing") = {1, 3, 4, 2};
//+
Physical Curve("Reflecting") = {5, 8, 7, 6};
//+
Physical Surface("Domain") = {9};
