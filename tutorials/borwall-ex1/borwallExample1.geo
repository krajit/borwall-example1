Point(newp) = {0,0,0};
Point(newp) = {0.3333,0,0};
Point(newp) = {0.6666,0,0};
Point(newp) = {1,0,0};
Point(newp) = {0.3333,1,0};
Point(newp) = {0.6666,1,0};
Point(newp) = {1,0.3333,0};
Point(newp) = {1,0.6666,0};
Point(newp) = {1,1,0};
Point(newp) = {0,1,0};
Point(newp) = {0,0.3333,0};
Point(newp) = {0,0.6666,0};


Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 7};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 9};
//+
Line(7) = {9, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {5, 10};
//+
Line(10) = {10, 12};
//+
Line(11) = {12, 11};
//+
Line(12) = {11, 1};
//+
Curve Loop(1) = {11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//+
Plane Surface(1) = {1};
Transfinite Curve{11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 33;
Transfinite Surface{1} = {1,4,9,10};
Recombine Surface{1};
//+
Extrude {0, 0, 1} {
  Surface{1}; 
  Layers{1};
  Recombine;
}
//+
Physical Surface("inlet") = {33,29,73};
Physical Surface("outlet") = {53};
Physical Surface("lowerWall") = {37,41,45,49};
Physical Surface("upperWall") = {69,65,61,57};//+
Physical Volume("body") = {1};
