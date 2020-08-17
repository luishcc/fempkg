cl__1 = 1;
Point(1) = {0, 0, 0, 0.1};
Point(2) = {5, 0, 0, 0.13};
Point(3) = {5, 1, 0, 0.09};
Point(4) = {0, 1, 0, 0.12};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};
Physical Line("dirichlet 0") = {1};
Physical Line("dirichlet 1") = {2};
Physical Line("dirichlet 2") = {3};
Physical Line("dirichlet 3") = {4};
Physical Surface(11) = {6};


