cl__1 = 1;
cl__2 = 0.3;
Point(3) = {0, 0.2, 0, 0.05};
Point(4) = {0, 0.3, 0, 0.05};
Point(7) = {1, 0.2, 0, 0.05};
Point(8) = {1, 0.3, 0, 0.05};


Line(3) = {3, 4};
Line(4) = {4, 8};
Line(8) = {7, 8};
Line(10) = {3, 7};


Line Loop(11) = {3, 4, -8, -10};
Plane Surface(12) = {11};


Physical Line("in hot") = {3};
Physical Line("out hot") = {8};
Physical Line("psi top") = {4};
Physical Line("psi bot") = {10};


Physical Surface(26) = {12};

