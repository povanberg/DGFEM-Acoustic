Point(1) = {0,0,0,0.5};
Point(2) = {1,0,0,1};
Point(3) = {1,1,0,1};
Point(4) = {0,1,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(10) = {1,2,3,4};
Plane Surface(20) = {10};
Physical Curve("dirichelet",100) = {1,2,3,4};
Physical Surface("domain",200) = {20};

Mesh 100;