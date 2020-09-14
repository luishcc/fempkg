
wall = 0.15;
hole = 0.05;

/* 
 *        4                 9L                 3
 *         o -------------------------------- o         
 *         |                                  |       
 *         |                       _    _     |       
 *      2L |                      / \    | L  |       Y         
 *         |                      \_/   _|    | 2L    ^
 *         |                                  |       |
 *         |                                  |       |
 *      1  o -------------------------------- o       o -----> X
 *         x(0,0)           9L                 2
 * */

D = 1.0;
xCenter=0.0;
yCenter=0.0;
radius=D/2.0;
pert=(00.0/100)*radius;
Point(1) = {xCenter, yCenter, 0,hole}; // center
Point(2) = {xCenter, yCenter+radius-pert, 0,hole}; // up
Point(3) = {xCenter, yCenter-radius+pert, 0,hole}; // down
Point(4) = {xCenter-radius-pert, yCenter, 0,hole}; // left
Point(5) = {xCenter+radius+pert, yCenter, 0,hole}; // right
Ellipse(1) = {2, 1, 1, 5};
Ellipse(2) = {5, 1, 1, 3};
Ellipse(3) = {3, 1, 1, 4};
Ellipse(4) = {4, 1, 1, 2};

Point(6)  = {-3.0*D, -3.0*D, 0.0,  wall}; // p1
Point(7)  = { 8.0*D, -3.0*D, 0.0,  wall}; // p2
Point(8)  = { 8.0*D,  3.0*D, 0.0,  wall}; // p3
Point(9)  = {-3.0*D,  3.0*D, 0.0,  wall}; // p4

Line(5) = {7, 6};
Line(6) = {8, 7};
Line(7) = {9, 8};
Line(8) = {6, 9};

Point(10) = {-1.0*D, -0.9*D, 0, hole};
Point(11) = { 1.4*D, -0.9*D, 0, hole};
Point(12) = { 1.4*D,  0.9*D, 0, hole};
Point(13) = {-1.0*D,  0.9*D, 0, hole};

Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};

//+
//Physical Line('wallNoSlip') = {5,7};
Physical Line('inlet') = {8};
Physical Line('outlet') = {6};
Physical Line('cylinder') = {3, 2, 1, 4}; 
//Physical Line('wallNoSlip') = {5}; 
//Physical Line('wallNoSlipConst') = {7}; 
Physical Line('bot') = {5}; 
Physical Line('top') = {7}; 
//Physical Line('interface') = {9,10,11,12}; 


//+
Line Loop(1) = {8, 7, 6, 5};
Line Loop(2) = {11, 12, 9, 10};
Plane Surface(1) = {1, 2};
Line Loop(3) = {1, 2, 3, 4};
Plane Surface(2) = {2, 3};
Physical Surface('inner') = {1};
Physical Surface('outer') = {2}; 



