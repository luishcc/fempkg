cl__1 = 1;
Point(1) = {0, 0, 0, 0.7};
Point(2) = {0, 5, 0, 0.7};
Point(3) = {7, 0, 0, 0.7};
Point(4) = {7, 5, 0, 0.7};
Point(5) = {3.5, 2.5, 0, 0.7};
Point(6) = {3.5, 1.5, 0, 0.3};
Point(7) = {2.633974596, 3, 0, 0.3};
Point(8) = {4.366025404, 3, 0, 0.3};


Point(9) = {0, 3, 0, 0.2};

Point(11) = {7, 3, 0, 0.2};



Line(1) = {1, 9};
Line(2) = {9,  2};
Line(4) = {2, 4};
Line(5) = {7, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(9) = {4, 11};
Line(11) = {11, 3};
Line(12) = {3, 1};
Line(8) = {9, 7};
Line(13) = {8, 11};


Line Loop(14) = {2, 4, 9, -13, 7, -8};
Plane Surface(15) = {14};
Line Loop(16) = {1, 8, 5, 6, 13, 11, 12};
Plane Surface(17) = {16};


Physical Line("neumann 0") = {1, 2, 4, 9, 11};
Physical Line("neumann conv") = {5, 7, 6};
Physical Line("neumann 1") = {12};
Physical Surface(21) = {15, 17};
