scale = 0.6;

Point(1) = {0, 0, 0, 1.5*scale};
Point(2) = {0, 10, 0, 1.5*scale};
Point(3) = {32.5, 10, 0, 1.5*scale};
Point(4) = {32.5, 5, 0, 1.0*scale};
Point(5) = {32.5, 0, 0, 1.5*scale};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Point(6) = {12, 4.5, 0, 0.1*scale};
Point(7) = {12, 5.5, 0, 0.1*scale};
Point(8) = {13, 5.5, 0, 0.1*scale};
Point(9) = {13, 4.5, 0, 0.1*scale};

Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 6};


Point(10) = {13, 5, 0, 0.1*scale};
Line(10) = {10, 4};


Line Loop(1) = {2, 3, 4, 5, 1};
Line Loop(2) = {6, 7, 8, 9};
Plane Surface(1) = {1, 2};

Field[1] = Attractor;
Field[1].EdgesList = {10};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.2*scale;
Field[2].LcMax = 0.8*scale;
Field[2].DistMin = 0.75;
Field[2].DistMax = 5;
Background Field = 2;

Physical Line("inlet") = {1};
Physical Line("outlet") = {3, 4};
Physical Line("top") = {2};
Physical Line("bot") = {5};
Physical Line("cylinder") = {6,7,8,9};

Physical Surface(6) = {1};
