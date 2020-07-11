cl__1 = 1;
Point(1) = {0, 0, 0, 0.7};
Point(2) = {0, 5, 0, 0.7};
Point(3) = {7, 0, 0, 0.7};
Point(4) = {7, 5, 0, 0.7};
Point(5) = {3.5, 2.65, 0, 0.7};
Point(6) = {4.936, 2.65, 0, 0.2};
Point(7) = {2.064, 2.65, 0, 0.2};
Point(8) = {3.5, 2.3, 0, 0.2};
Point(13) = {3.5, 3, 0, 0.2};

Point(9) = {0, 3, 0, 0.2};
Point(11) = {7, 3, 0, 0.2};
Point(12) = {3.5, 3.01, 0, 0.4};

Line(1) = {1, 9};
Line(2) = {9, 2};
Line(4) = {2, 4};
Line(5) = {9, 12};

Ellipse(31) = {13, 5, 7, 7};
Ellipse(32) = {7, 5, 6, 8};
Ellipse(33) = {8, 5, 6, 6};
Ellipse(34) = {6, 5, 7, 13};

Line(6) = {12, 11};
Line(9) = {4, 11};
Line(11) = {11, 3};
Line(12) = {3, 1};



Line Loop(42) = {2, 4, 9, -6, -5};
Plane Surface(43) = {42};
Line Loop(44) = {5, 6, 11, 12, 1};
Line Loop(45) = {32, 33, 34, 31};
Plane Surface(46) = {44, 45};


Physical Line("neumann 0") = {1, 2, 4, 9, 11};
Physical Line("neumann conv") = {33, 32, 31, 34};
Physical Line("neumann 1") = {12};

Physical Surface(47) = {43, 46};
