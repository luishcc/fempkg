
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 0.5, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {2, 1, 0, 1.0};
Point(5) = {2, 0.5, 0, 1.0};
Point(6) = {2, 0, 0, 1.0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};


Line Loop(1) = {2, 3, 4, -7};
Plane Surface(1) = {1};
Line Loop(2) = {1, 7, 5, 6};
Plane Surface(2) = {2};


Field[3] = Attractor;
Field[3].NNodesByEdge = 100;
Field[3].EdgesList = {7};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.03;
Field[4].LcMax = 0.3;
Field[4].DistMin = 0.02;
Field[4].DistMax = 0.9;
Background Field = 4;



Physical Line("in") = {2};
Physical Line("out") = {4};
Physical Line("top") = {3};
Physical Line("bot") = {7};
Physical Line("wallbot") = {6};

Physical Surface(101) = {1, 2};
