cl__1 = 1;

Point(1) = {0, 0, 0, 0.06};
Point(2) = {5, 0, 0, 0.06};
Point(3) = {5, 1, 0, 0.06};
Point(4) = {0, 1, 0, 0.06};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

Line Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};

Physical Line("inflow") = {1};
Physical Line("top") = {2};
Physical Line("outflow") = {3};
Physical Line("bot") = {4};

Physical Surface(11) = {6};


