cl__1 = 0.2;
Point(1) = {0, 0, 0, 0.08};
Point(2) = {1, 0, 0, 0.08};
Point(3) = {1, 1, 0, 0.08};
Point(4) = {0, 1, 0, 0.08};
Point(5) = {0, 0.5, 0, 0.08}; 
Point(6) = {1, 0.5, 0, 0.08};

Line(1) = {1, 5};
Line(2) = {5, 4};
Line(3) = {4, 3};
Line(4) = {3, 6};
Line(5) = {6, 2};
Line(6) = {2, 1};
Line(7) = {5, 6};

Line Loop(11) = {1, 7, 5, 6};
Plane Surface(21) = {11};

Line Loop(12) = {2, 3, 4, -7};
Plane Surface(22) = {12};

Physical Line("dirichlet 0") = {6};
Physical Line("dirichlet 1") = {3};
Physical Line(100) = {1, 2, 4, 5, 7};

Physical Surface(101) = {21};
Physical Surface(102) = {22};
