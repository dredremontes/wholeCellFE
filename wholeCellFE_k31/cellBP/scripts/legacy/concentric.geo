// Gmsh project created on Tue May 30 11:17:28 2023
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 9.4, 0, 2*Pi};
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {2};
//+
Plane Surface(1) = {3, 4};
//+
Physical Curve(5) = {2};
//+
Physical Curve(6) = {1};
//+
Curve Loop(5) = {2};
//+
Curve Loop(6) = {2};
//+
Surface(2) = {5, 6};
//+
Curve Loop(7) = {2};
//+
Plane Surface(2) = {7};

//+
Field[1] = Ball;
//+
Field[1].Radius = 5;
//+
Field[1].Thickness = 5;
//+
Field[1].VIn = 1;
//+
Field[1].VOut = 1.5;
//+
Background Field = 1;
//+
Field[1].VIn = 0.8;
//+
Field[1].VOut = 1;
//+
Field[1].VIn = 0.7;
//+
Field[1].Radius = 7;
//+
Field[1].Thickness = 3;
//+
Field[1].Radius = 3;
//+
Field[1].Thickness = 7;
//+
Field[1].VIn = 1;
//+
Field[1].VOut = 0.7;
//+
Field[1].VOut = 0.75;
//+
Field[1].VOut = 0.8;
//+
Field[1].VOut = 0.9;
//+
Field[1].VIn = 2;
//+
Field[1].VOut = 1.5;
//+
Field[1].VOut = 1.3;
//+
Field[1].Radius = 5;
//+
Field[1].Thickness = 5;
//+
Field[1].Radius = 1;
//+
Field[1].Thickness = 9;
//+
Field[1].VOut = 1.5;
//+
Field[1].VIn = 1.5;
//+
Field[1].VOut = 1;
