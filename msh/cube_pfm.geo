/*********************************************************************
 *
 *  Mode I fracture - Opening mode
 *
 *********************************************************************/

L = 10;
L_crack = 1e-3;

lc = 1;
lc_crack = 1e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc_crack};
Point(3) = {0, L, 0, lc};
Point(4) = {L, L, 0, lc_crack};
Point(5) = {0,   (L-L_crack)/2, 0, lc};
Point(6) = {0,   (L+L_crack)/2, 0, lc};
Point(7) = {L/2, L/2, 0, lc_crack};


Line(1)  = {1,2};
Line(2)  = {2,4};
Line(3)  = {4,3};
Line(4)  = {3,6};
Line(5)  = {6,7};
Line(6)  = {7,5};
Line(7)  = {5,1};

Line Loop(8) = {  1,  2,  3,  4,  5,  6,  7};

Plane Surface(1) = {8};

Physical Line("disp_fix",2) = {1};
Physical Line("disp",3) = {3};
Physical Surface("surface",1) = {1};

Mesh.MshFileVersion = 2;
Mesh 2;
