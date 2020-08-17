cl__1 = 1;
cl__2 = 0.3;
Point(1) = {0, 0, 0, 0.05};
Point(2) = {0, 0.1, 0, 0.05};
Point(3) = {0, 0.2, 0, 0.05};
Point(4) = {0, 0.3, 0, 0.05};
Point(5) = {1, 0, 0, 0.05};
Point(6) = {1, 0.1, 0, 0.05};
Point(7) = {1, 0.2, 0, 0.05};
Point(8) = {1, 0.3, 0, 0.05};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 8};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {2, 6};
Line(10) = {3, 7};


Line Loop(11) = {3, 4, -8, -10};
Plane Surface(12) = {11};
Line Loop(13) = {1, 9, -6, -5};
Plane Surface(14) = {13};
Line Loop(15) = {10, -7, -9, 2};
Plane Surface(16) = {15};


Physical Line("in hot") = {3};
Physical Line("in cold") = {6};
Physical Line("out hot") = {8};
Physical Line("out cold") = {1};
Physical Line(22) = {4};
Physical Line(23) = {10};
Physical Line(24) = {9};
Physical Line(25) = {5};

Physical Surface(26) = {12};
Physical Surface(27) = {14};
Physical Surface(28) = {16};

