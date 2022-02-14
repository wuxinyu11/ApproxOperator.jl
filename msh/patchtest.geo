
a = 1.0;
b = 1.0;
n = 11;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  a,   b, 0.0};
Point(4) = {0.0,   b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4} = n;

// Physical Curve("Γᵗ₁") = {2};
// Physical Curve("Γᵗ₂") = {3};
Physical Curve("Γᵍ") = {1,2,3,4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 9;
Mesh.MshFileVersion = 2;
Mesh 2;
// RecombineMesh;
