cl__1 = 1;
cl__2 = 0.3;
Point(1) = {0, 0, 0, 0.1};
Point(2) = {2, 0, 0, 0.1};
Point(3) = {2, 0.25, 0, 0.03};
Point(4) = {0, 0.25, 0, 0.03};
Point(5) = {0, 0.75, 0, 0.1};
Point(6) = {2, 0.75, 0, 0.1};

Line(1) = {1, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {1, 2};
Line(5) = {4, 3};
Line(6) = {6, 3};
Line(7) = {3, 2};

Line Loop(1) = {2, 3, 6, -5};
Plane Surface(1) = {1};
Line Loop(2) = {1, 5, 7, -4};
Plane Surface(2) = {2};

Physical Line("neumann 0") = {1, 2};
Physical Line("neumann 1") = {7, 6};
Physical Line("dirichlet 1") = {4};
Physical Line("dirichlet 0") = {3};


Physical Surface(5) = {1};
Physical Surface(6) = {2};
