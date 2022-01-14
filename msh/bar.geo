
a = 1.0;
n = 10;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {1.0, 0.0, 0.0};

Line(1) = {1,2};

Transfinite Curve{1} = n;

Physical Point("EBC") = {1};
Physical Point("NBC") = {2};
Physical Curve("Domain") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;
