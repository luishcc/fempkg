cl__1 = 1;
cl__2 = 0.3;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {0, 1.5, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};
Point(5) = {5, 2, 0, 1.0};
Point(6) = {5, 1.5, 0, 1.0};
Point(7) = {5, 1, 0, 1.0};
Point(8) = {5, 0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line(9) = {2, 7};
Line(10) = {3, 6};

Line Loop(1) = {3, 4, 5, -10};
Plane Surface(1) = {1};
Line Loop(2) = {2, 10, 6, -9};
Plane Surface(2) = {2};
Line Loop(3) = {8, 1, 9, 7};
Plane Surface(3) = {3};

Field[3] = Attractor;
Field[3].NNodesByEdge = 100;
Field[3].EdgesList = {9, 4};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.04;
Field[4].LcMax = 0.2;
Field[4].DistMin = 0.02;
Field[4].DistMax = 0.9;
Background Field = 4;


Physical Line("top fluid") = {4};
Physical Line("interface") = {9};
Physical Line("bot solid") = {8};
Physical Line("left solid") = {1};
Physical Line("right solid") = {7};
Physical Line("in") = {3, 2};
Physical Line("out") = {5, 6};

Physical Surface(8) = {1, 2, 3};



