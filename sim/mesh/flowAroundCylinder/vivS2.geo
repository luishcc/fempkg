
Point(1) = {5, 0, 0, 1.5};
Point(2) = {5, 10, 0, 1.5};
Point(3) = {25, 10, 0, 1.5};
Point(4) = {25, 5, 0, 0.2};
Point(5) = {25, 0, 0, 1.5};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Point(6) = {12.5, 5, 0, 0.1};
Point(7) = {13, 5, 0, 0.1};

Line(10) = {7, 4};


Circle(6) = {7, 6, 7};

Line Loop(1) = {2, 3, 4, 5, 1};
Line Loop(2) = {6};
Plane Surface(1) = {1, 2};

Field[1] = Attractor;
Field[1].EdgesList = {10};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.1;
Field[2].LcMax = 0.8;
Field[2].DistMin = 0.5;
Field[2].DistMax = 5;
Background Field = 2;

Physical Line("inlet") = {1};
Physical Line("outlet") = {3, 4};
Physical Line("top") = {2};
Physical Line("bot") = {5};
Physical Line("cylinder") = {6};

Physical Surface(6) = {1};
