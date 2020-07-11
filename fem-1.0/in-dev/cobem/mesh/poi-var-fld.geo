cl__1 = 1;
cl__2 = 0.3;

Point(1) = {0, 0, 0, 0.5};
Point(2) = {0, 1, 0, 0.5};
Point(3) = {5, 1, 0, 0.3};
Point(4) = {5, 0, 0, 0.3};
Point(5) = {0, 0.5, 0, 0.5};
Point(6) = {5, 0.5, 0, 0.5};

Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};

Line Loop(11) = {1, 7, 5, 6};
Plane Surface(12) = {11};
Line Loop(21) = {2, 3, 4, -7};
Plane Surface(22) = {21};


Field[3] = Attractor;
Field[3].NNodesByEdge = 100;
Field[3].EdgesList = {3, 6};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.03;
Field[4].LcMax = 0.3;
Field[4].DistMin = 0.02;
Field[4].DistMax = 0.9;
Background Field = 4;


Physical Line("in") = {1, 2};
Physical Line("out") = {4, 5};
Physical Line("top") = {3};
Physical Line("bot") = {6};

Physical Surface(26) = {12};
Physical Surface(27) = {22};


