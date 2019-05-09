//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("Absorbing") = {1};
//+
Physical Volume("Domain") = {1};
