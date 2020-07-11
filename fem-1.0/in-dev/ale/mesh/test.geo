cl__1 = 1;
cl__2 = 0.3;

Point(1) = {0, 0, 0, 2.8};
Point(2) = {0, 1, 0, 0.5};
Point(3) = {1, 1, 0, 1.3};
Point(4) = {1, 0, 0, 0.1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};

Physical Line("dirichlet 0") = {1, 3, 2, 4};

Physical Surface(26) = {12};

