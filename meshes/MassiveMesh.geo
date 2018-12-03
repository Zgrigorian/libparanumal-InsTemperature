//+
Point(1) = {-10, -10, 0, 1.0};
//+
Point(2) = {-10, 10, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {10, -10, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Wall") = {1, 3};
//+
Physical Curve("Inlet") = {4};
//+
Physical Curve("outlet") = {2};
