// Define the radius of the circle
r = 10;

// Define the geometry of the circle
Point(1) = {0, 0, 0, 1.0};
Point(2) = {r, 0, 0, 1.0};
Point(3) = {0, r, 0, 1.0};
Point(4) = {-r, 0, 0, 1};
Point(5) = {0, -r, 0, 1};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1) = {2,3,4,1};
Plane Surface(1) = {1};