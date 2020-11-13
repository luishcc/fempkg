cl__1 = 1;
cl__2 = 0.;
h=0.16;

Point(1) = {0, 0, 0, 0.4*h};
Point(2) = {0, 4, 0, 0.4*h};
Point(3) = {4, 4, 0, 0.4*h};
Point(4) = {4, 0, 0, 0.4*h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};

Physical Line("in") = {1};
Physical Line("out") = {3};
Physical Line("top") = {2};
Physical Line("bot") = {4};

Physical Surface(26) = {12};

