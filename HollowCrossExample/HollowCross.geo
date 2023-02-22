// Gmsh project created on Thu Sep 22 11:56:23 2022
// Hollow cross Fig.9(a)

//---------------------------------------------------------------------------
// All the nodes
h1 = 1.0;
Point(1) = {0,   0,    0,   h1};
Point(2) = {11,  0,    0,   h1};
Point(3) = {0,   6,    0,   h1};
Point(4) = {5,   6,    0,   h1};
Point(5) = {5,   11,   0,   h1};
Point(6) = {11,  11,   0,   h1};
Point(7) = {11,  16,   0,   h1};
Point(8) = {5,   22,   0,   h1};
Point(9) = {16,  16,   0,   h1};
Point(10) = {16,  22,   0,   h1};
Point(11) = {27,  16,   0,   h1};
Point(12) = {21,  22,   0,   h1};
Point(13) = {27,  27,   0,   h1};
Point(14) = {21,  27,   0,   h1};
Point(15) = {27,  38,   0,   h1};
Point(16) = {21,  32,   0,   h1};
Point(17) = {16,  38,   0,   h1};
Point(18) = {16,  32,   0,   h1};
Point(19) = {11,  38,   0,   h1};
Point(20) = {5,   32,   0,   h1};
Point(21) = {11,  43,   0,   h1};
Point(22) = {5,   43,   0,   h1};
Point(23) = {11,  54,   0,   h1};
Point(24) = {5,   48,   0,   h1};
Point(25) = {0,   54,   0,   h1};
Point(26) = {0,   48,   0,   h1};

Point(27) = {-11,  0,    0,   h1};
Point(28) = {-5,   6,    0,   h1};
Point(29) = {-5,   11,   0,   h1};
Point(30) = {-11,  11,   0,   h1};
Point(31) = {-11,  16,   0,   h1};
Point(32) = {-5,   22,   0,   h1};
Point(33) = {-16,  16,   0,   h1};
Point(34) = {-16,  22,   0,   h1};
Point(35) = {-27,  16,   0,   h1};
Point(36) = {-21,  22,   0,   h1};
Point(37) = {-27,  27,   0,   h1};
Point(38) = {-21,  27,   0,   h1};
Point(39) = {-27,  38,   0,   h1};
Point(40) = {-21,  32,   0,   h1};
Point(41) = {-16,  38,   0,   h1};
Point(42) = {-16,  32,   0,   h1};
Point(43) = {-11,  38,   0,   h1};
Point(44) = {-5,   32,   0,   h1};
Point(45) = {-11,  43,   0,   h1};
Point(46) = {-5,   43,   0,   h1};
Point(47) = {-11,  54,   0,   h1};
Point(48) = {-5,   48,   0,   h1};

//---------------------------------------------------------------------------
// All the lines
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Line(5) = {5, 4};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 2};
//+
Line(8) = {7, 6};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 5};
//+
Line(11) = {7, 9};
//+
Line(12) = {9, 10};
//+
Line(13) = {10, 8};
//+
Line(14) = {10, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 9};
//+
Line(17) = {11, 13};
//+
Line(18) = {13, 14};
//+
Line(19) = {14, 12};
//+
Line(20) = {14, 16};
//+
Line(21) = {16, 15};
//+
Line(22) = {15, 13};
//+
Line(23) = {15, 17};
//+
Line(24) = {17, 18};
//+
Line(25) = {18, 16};
//+
Line(26) = {18, 20};
//+
Line(27) = {20, 19};
//+
Line(28) = {19, 17};
//+
Line(29) = {21, 19};
//+
Line(30) = {21, 22};
//+
Line(31) = {22, 20};
//+
Line(32) = {22, 24};
//+
Line(33) = {24, 23};
//+
Line(34) = {23, 21};
//+
Line(35) = {23, 25};
//+
Line(36) = {25, 26};
//+
Line(37) = {26, 24};
//+
Line(38) = {25, 47};
//+
Line(39) = {47, 48};
//+
Line(40) = {48, 26};
//+
Line(41) = {48, 46};
//+
Line(42) = {46, 45};
//+
Line(43) = {45, 47};
//+
Line(44) = {46, 44};
//+
Line(45) = {44, 43};
//+
Line(46) = {43, 45};
//+
Line(47) = {43, 41};
//+
Line(48) = {41, 42};
//+
Line(49) = {42, 44};
//+
Line(50) = {42, 40};
//+
Line(51) = {40, 39};
//+
Line(52) = {39, 41};
//+
Line(53) = {39, 37};
//+
Line(54) = {37, 38};
//+
Line(55) = {38, 40};
//+
Line(56) = {37, 35};
//+
Line(57) = {35, 36};
//+
Line(58) = {36, 38};
//+
Line(59) = {35, 33};
//+
Line(60) = {33, 34};
//+
Line(61) = {34, 36};
//+
Line(62) = {34, 32};
//+
Line(63) = {32, 29};
//+
Line(64) = {29, 30};
//+
Line(65) = {30, 31};
//+
Line(66) = {31, 33};
//+
Line(67) = {31, 32};
//+
Line(68) = {30, 27};
//+
Line(69) = {27, 1};
//+
Line(70) = {29, 28};
//+
Line(71) = {28, 3};
//+
Line(72) = {28, 27};

//---------------------------------------------------------------------------
// This section describes the "Plane Surfaces", i.e., the 2D surfaces for meshing

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, -2, -7, -6};
Plane Surface(2) = {2};

Curve Loop(3) = {6, -8, 9, 10};
Plane Surface(3) = {3};

Curve Loop(4) = {11, 12, 13, -9};
Plane Surface(4) = {4};

Curve Loop(5) = {-15, -14, -12, -16};
Plane Surface(5) = {5};

Curve Loop(6) = {15, 17, 18, 19};
Plane Surface(6) = {6};

Curve Loop(7) = {-18, -22, -21, -20};
Plane Surface(7) = {7};

Curve Loop(8) = {21, 23, 24, 25};
Plane Surface(8) = {8};

Curve Loop(9) = {-24, -28, -27, -26};
Plane Surface(9) = {9};

Curve Loop(10) = {27, -29, 30, 31};
Plane Surface(10) = {10};

Curve Loop(11) = {-30, -34, -33, -32};
Plane Surface(11) = {11};

Curve Loop(12) = {37, 33, 35, 36};
Plane Surface(12) = {12};

Curve Loop(13) = {-36, 38, 39, 40};
Plane Surface(13) = {13};

Curve Loop(14) = {-39, -43, -42, -41};
Plane Surface(14) = {14};

Curve Loop(15) = {42, -46, -45, -44};
Plane Surface(15) = {15};

Curve Loop(16) = {45, 47, 48, 49};
Plane Surface(16) = {16};

Curve Loop(17) = {-48, -52, -51, -50};
Plane Surface(17) = {17};

Curve Loop(18) = {51, 53, 54, 55};
Plane Surface(18) = {18};

Curve Loop(19) = {-54, 56, 57, 58};
Plane Surface(19) = {19};

Curve Loop(20) = {-57, 59, 60, 61};
Plane Surface(20) = {20};

Curve Loop(21) = {-60, -66, 67, -62};
Plane Surface(21) = {21};

Curve Loop(22) = {-67, -65, -64, -63};
Plane Surface(22) = {22};

Curve Loop(23) = {64, 68, -72, -70};
Plane Surface(23) = {23};

Curve Loop(24) = {72, 69, -4, -71};
Plane Surface(24) = {24};

//---------------------------------------------------------------------------
// Creates physical surfaces. The number in () represents the material id.
// The expression list on the right hand side is the list
// of elementary surfaces created above.  This is what makes our mesh 2D.
Physical Surface(1) = {1, 5, 8,  12, 16, 21};
Physical Surface(2) = {4, 9, 13, 17, 20, 24};
Physical Surface(3) = {2, 6, 10, 15, 19, 23};
Physical Surface(4) = {3, 7, 11, 14, 18, 22};

//---------------------------------------------------------------------------
// Parameters for the meshing
Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 2.0;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
