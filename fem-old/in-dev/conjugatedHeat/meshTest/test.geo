cl__1 = 1;
cl__2 = 0.3;
Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};
Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0, 0.5};
Point(5) = {0, 2, 0, 0.5};
Point(6) = {1, 2, 0, 0.5};

Line(1) = {1, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {1, 2};
Line(5) = {4, 3};
Line(6) = {6, 3};
Line(7) = {3, 2};

Line Loop(9) = {2, 3, 6, -5};
Plane Surface(9) = {9};
Line Loop(11) = {4, -7, -5, -1};
Plane Surface(11) = {11};



Physical Line("dirichlet 0") = {3};
Physical Line(31) = {5};
Physical Line("neumann 0") = {1, 2};
Physical Line("neumann 0") = {7, 6};
Physical Line("dirichlet 1") = {4};

Physical Surface(21) = {9};
Physical Surface(22) = {11};
