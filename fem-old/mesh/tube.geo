cl__1 = 0.1;
Point(1) = {0.5, 0.5, 0, 0.1};
Point(2) = {0.5, 0, 0, 0.1};
Point(3) = {0.5, 0.3, 0, 0.1};
Circle(1) = {2, 1, 2};
Circle(2) = {3, 1, 3};
Line Loop(5) = {1, -2};

Physical Line("Out") = {1};
Physical Line("In") = {2};

Line Loop(6) = {2};
Plane Surface(7) = {5, 6};
Physical Surface(8) = {7};
