cl__1 = 0.1;
Point(1) = {0, 0, 0, 0.1};
Point(2) = {4, 0, 0, 0.1};
Point(3) = {4, 1, 0, 0.1};
Point(4) = {4, 2, 0, 0.3};
Point(5) = {0, 2, 0, 0.3};
Point(6) = {0, 1, 0, 0.1};

Line(1) = {1, 6};
Line(2) = {6, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};
Line(6) = {2, 1};
Line(7) = {6, 3};

Line Loop(9) = {2, 3, 4, -7};
Plane Surface(9) = {9};
Line Loop(11) = {1, 7, 5, 6};
Plane Surface(11) = {11};

Physical Line("dirichlet 0") = {6};
Physical Line("dirichlet 1") = {3};
Physical Line(14) = {2, 1, 5, 4};

Physical Surface(13) = {9, 11};

