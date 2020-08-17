
Point(1) = {0, 0, 0, 1.2};
Point(2) = {0, 10, 0, 1.2};
Point(3) = {32.5, 10, 0, 1.2};
Point(4) = {32.5, 0, 0, 1.2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {12.5, 5, 0, 0.1};
Point(6) = {13, 5, 0, 0.1};

Circle(5) = {6, 5, 6};

Point(7) = {12.5, 5, 0, 0.8};
Point(8) = {15.5, 5, 0, 0.8};

Circle(6) = {8, 7, 8};

Line Loop(1) = {2, 3, 4, 1};
Line Loop(2) = {6};
Plane Surface(1) = {1, 2};

Line Loop(3) = {6};
Line Loop(4) = {5};
Plane Surface(2) = {3, 4};

Physical Line("inlet") = {1};
Physical Line("outlet") = {3};
Physical Line("top") = {2};
Physical Line("bot") = {4};
Physical Line("cylinder") = {5};
Physical Line("mesh") = {6};

Physical Surface(6) = {1,2};
