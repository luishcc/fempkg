

Point(2) = {0, 0.5, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {2, 1, 0, 1.0};
Point(5) = {2, 0.5, 0, 1.0};


Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(7) = {2, 5};


Line Loop(1) = {2, 3, 4, -7};
Plane Surface(1) = {1};

Field[3] = Attractor;
Field[3].NNodesByEdge = 100;
Field[3].EdgesList = {7};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.03;
Field[4].LcMax = 0.3;
Field[4].DistMin = 0.02;
Field[4].DistMax = 0.9;
Background Field = 4;


Physical Line("neumann 0") = {2};
Physical Line("neumann") = {4};
Physical Line("dirichlet 1") = {3};
Physical Line("dirichlet 0") = {7};

Physical Surface(101) = {1};
