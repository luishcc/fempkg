cl__1 = 1;
cl__2 = 0.3;
Point(1) = {0, 0, 0, 0.05};
Point(2) = {0, 0.1, 0, 0.05};
Point(5) = {1, 0, 0, 0.05};
Point(6) = {1, 0.1, 0, 0.05};

Line(1) = {1, 2};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(9) = {2, 6};


Line Loop(13) = {1, 9, -6, -5};
Plane Surface(14) = {13};


Physical Line("in cold") = {6};
Physical Line("out cold") = {1};
Physical Line("psi top") = {9};
Physical Line("psi bot") = {5};


Physical Surface(27) = {14};


