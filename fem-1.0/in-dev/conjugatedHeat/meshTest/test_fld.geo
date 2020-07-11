cl__1 = 1;
cl__2 = 0.3;

Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0, 0.5};
Point(5) = {0, 2, 0, 0.5};
Point(6) = {1, 2, 0, 0.5};


Line(2) = {4, 5};
Line(3) = {5, 6};
Line(5) = {4, 3};
Line(6) = {6, 3};

Line Loop(9) = {2, 3, 6, -5};
Plane Surface(9) = {9};

Physical Line("dirichlet 0") = {3};
Physical Line("neumann 0") = {2};
Physical Line("neumann 1") = {6};
Physical Line("dirichlet 1") = {5};

Physical Surface(21) = {9};

