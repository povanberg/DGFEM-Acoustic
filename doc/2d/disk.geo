//+
SetFactory("OpenCASCADE");
Disk(1) = {-0, 0, 0, 6, 6};
//+
Physical Surface("Domain") = {1};
//+
Physical Curve("Reflecting") = {1};
