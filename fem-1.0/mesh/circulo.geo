cl__1 = 1;
Point(1) = {0, 0, 0, 0.7};
Point(2) = {0, 5, 0, 0.7};
Point(3) = {7, 0, 0, 0.7};
Point(4) = {7, 5, 0, 0.7};
Point(5) = {3.5, 2.5, 0, 0.7};
Point(6) = {3.5, 2, 0, 0.2};
Point(7) = {3.5, 3, 0, 0.2};



Point(9) = {0, 3.01, 0, 0.2};
Point(11) = {7, 3.01, 0, 0.2};
Point(12) = {3.5, 3.01, 0, 0.2};


Line(1) = {1, 9};
Line(2) = {9, 2};
Line(4) = {2, 4};
Line(5) = {9, 12};
Circle(28) = {7, 5, 7};
Line(6) = {12, 11};
Line(9) = {4, 11};
Line(11) = {11, 3};
Line(12) = {3, 1};

Line Loop(29) = {2, 4, 9, -6, -5};
Plane Surface(30) = {29};
Line Loop(31) = {5, 6, 11, 12, 1};
Line Loop(32) = {28};
Plane Surface(33) = {31, 32};


Physical Line("neumann 0") = {1, 2, 4, 9, 11};
Physical Line("neumann conv") = {28};
Physical Line("neumann 1") = {12};

Physical Surface(37) = {30, 33};
