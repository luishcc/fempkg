cl__1 = 0.2;
Point(1) = {0, 0, 0, 0.08};
Point(2) = {3, 0, 0, 0.08};
Point(3) = {3, 2, 0, 0.08};
Point(4) = {0, 2, 0, 0.08};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Line("dirichlet 0") = {4, 2, 1};
Physical Line("dirichlet 1") = {3};
Physical Surface(9) = {6};
