cl__1 = 1;
Point(2) = {0, 0, 0, 0.1};
Point(3) = {0, 1, 0, 0.1};
Point(4) = {3, 1, 0, 0.1};
Point(5) = {3, 0, 0, 0.1};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Line Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};

Physical Line("dirichlet 500") = {2};
Physical Line("neumann 0") = {3};
Physical Line("neumann 1") = {4};
Physical Line("dirichlet 100") = {1};

Physical Surface(11) = {6};
