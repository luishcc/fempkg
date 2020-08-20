
Point(1) = {0, 0, 0, 3};
Point(2) = {0, 5, 0, 3};
Point(3) = {5, 5, 0, 3};
Point(4) = {5, 0, 0, 3};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {2.5, 2.5, 0, 1};
Point(6) = {3, 2.5, 0, 1};

Circle(5) = {6, 5, 6};

Line Loop(1) = {2, 3, 4, 1};
Line Loop(2) = {5};
Plane Surface(1) = {1, 2};

Physical Line("inlet") = {1};
Physical Line("outlet") = {3};
Physical Line("top") = {2};
Physical Line("bot") = {4};
Physical Line("cylinder") = {5};

Physical Surface("inner") = {1};
