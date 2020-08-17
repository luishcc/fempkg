cl__1 = 1;
Point(1) = {0, 0, 0, 0.70};
Point(2) = {0, 5.00, 0, 0.70};
Point(3) = {7.00, 0, 0, 0.70};
Point(4) = {7.00, 5.00, 0, 0.70};
Point(5) = {4.00, 2.00, 0, 0.30};
Point(6) = {3.00, 2.00, 0, 0.30};
Point(7) = {3.00, 3.000, 0, 0.30};
Point(8) = {4.00, 3.00, 0, 0.30};


Point(9) = {0, 3.00, 0, 0.20};
Point(11) = {7.00, 3.00, 0, 0.20};



Line(1) = {1, 9};
Line(2) = {9, 2};
Line(4) = {2, 4};
Line(5) = {7, 6};
Line(6) = {5, 8};
Line(7) = {8, 7};
Line(9) = {4, 11};
Line(11) = {11, 3};
Line(12) = {3, 1};
Line(8) = {9, 7};
Line(13) = {8, 11};

Line(15) = {5, 6};


Line Loop(16) = {2, 4, 9, -13, 7, -8};
Plane Surface(17) = {16};
Line Loop(18) = {1, 8, 5, -15, 6, 13, 11, 12};
Plane Surface(19) = {18};

Physical Line("neumann 0") = {1, 2, 4, 9, 11};
Physical Line("neumann cc") = {5, 15, 6, 7};
Physical Line("neumann 1") = {12};


Physical Surface(20) = {17, 19};
