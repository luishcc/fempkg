// Gmsh project created on Tue Apr 16 14:35:56 2019
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 1.0};
//+
Point(3) = {5, 1, 0, 1.0};
//+
Point(4) = {5, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Line(1) = {5, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 5};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Line("neumann 0") = {1, 3};
//+
Physical Line("dirichlet 1") = {4};
//+
Physical Line("dirichlet 0") = {2};
//+
Physical Surface(4) = {1};
