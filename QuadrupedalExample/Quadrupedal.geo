// Gmsh project created on Thu Oct 5 13:09:23 2022
// quadrupedal Fig.2(m) in "Printing ferromagnetic domains for untethered
// fast-transforming soft materials"

h1 = 1.0;
//+
Point(1) = {-9, -3, 0, h1};
//+
Point(2) = {-3, -3, 0, h1};
//+
Point(3) = {3, -3, 0, h1};
//+
Point(4) = {9, -3, 0, h1};
//+
Point(5) = {9, 3, 0, h1};
//+
Point(6) = {3, 3, 0, h1};
//+
Point(7) = {-3, 3, 0, h1};
//+
Point(8) = {-9, 3, 0, h1};
//+
Point(9) = {3, 9, 0, h1};
//+
Point(10) = {-3, 9, 0, h1};
//+
Point(11) = {-3, -9, 0, h1};
//+
Point(12) = {3, -9, 0, h1};
//+
Point(13) = {15, -3, 0, h1};
//+
Point(14) = {21, -3, 0, h1};
//+
Point(15) = {27, -3, 0, h1};
//+
Point(16) = {27, 3, 0, h1};
//+
Point(17) = {21, 3, 0, h1};
//+
Point(18) = {15, 3, 0, h1};
//+
Point(19) = {21, 9, 0, h1};
//+
Point(20) = {15, 9, 0, h1};
//+
Point(21) = {15, -9, 0, h1};
//+
Point(22) = {21, -9, 0, h1};
//+
Point(23) = {-27, -3, 0, h1};
//+
Point(24) = {-21, -3, 0, h1};
//+
Point(25) = {-15, -3, 0, h1};
//+
Point(26) = {-15, 3, 0, h1};
//+
Point(27) = {-21, 3, 0, h1};
//+
Point(28) = {-27, 3, 0, h1};
//+
Point(29) = {-15, 9, 0, h1};
//+
Point(30) = {-21, 9, 0, h1};
//+
Point(31) = {-21, -9, 0, h1};
//+
Point(32) = {-15, -9, 0, h1};
//+
Point(33) = {-9, 15, 0, h1};
//+
Point(34) = {-3, 15, 0, h1};
//+
Point(35) = {3, 15, 0, h1};
//+
Point(36) = {9, 15, 0, h1};
//+
Point(37) = {9, 21, 0, h1};
//+
Point(38) = {3, 21, 0, h1};
//+
Point(39) = {-3, 21, 0, h1};
//+
Point(40) = {-9, 21, 0, h1};
//+
Point(41) = {3, 27, 0, h1};
//+
Point(42) = {-3, 27, 0, h1};
//+
Point(43) = {-9, -21, 0, h1};
//+
Point(44) = {-3, -21, 0, h1};
//+
Point(45) = {3, -21, 0, h1};
//+
Point(46) = {9, -21, 0, h1};
//+
Point(47) = {9, -15, 0, h1};
//+
Point(48) = {3, -15, 0, h1};
//+
Point(49) = {-3, -15, 0, h1};
//+
Point(50) = {-9, -15, 0, h1};
//+
Point(51) = {-3, -27, 0, h1};
//+
Point(52) = {3, -27, 0, h1};
//+
Point(53) = {0, 0, 0, h1};
//+
Point(54) = {3, 0, 0, h1};
//+
Point(55) = {0, 3, 0, h1};
//+
Point(56) = {-3, 0, 0, h1};
//+
Point(57) = {0, -3, 0, h1};

//+
Line(1) = {23, 24};
//+
Line(2) = {24, 25};
//+
Line(3) = {25, 1};
//+
Line(4) = {1, 2};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 13};
//+
Line(8) = {13, 14};
//+
Line(9) = {14, 15};
//+
Line(10) = {28, 27};
//+
Line(11) = {27, 26};
//+
Line(12) = {26, 8};
//+
Line(13) = {8, 7};
//+
Line(15) = {6, 5};
//+
Line(16) = {5, 18};
//+
Line(17) = {18, 17};
//+
Line(18) = {17, 16};
//+
Line(19) = {23, 28};
//+
Line(20) = {24, 27};
//+
Line(21) = {25, 26};
//+
Line(22) = {1, 8};
//+
Line(25) = {4, 5};
//+
Line(26) = {13, 18};
//+
Line(27) = {14, 17};
//+
Line(28) = {15, 16};
//+
Line(29) = {30, 29};
//+
Line(30) = {29, 26};
//+
Line(31) = {30, 27};
//+
Line(32) = {24, 31};
//+
Line(33) = {31, 32};
//+
Line(34) = {32, 25};
//+
Line(35) = {20, 19};
//+
Line(36) = {19, 17};
//+
Line(37) = {20, 18};
//+
Line(38) = {13, 21};
//+
Line(39) = {21, 22};
//+
Line(40) = {22, 14};
//+
Line(41) = {2, 11};
//+
Line(42) = {11, 12};
//+
Line(43) = {12, 3};
//+
Line(44) = {10, 7};
//+
Line(45) = {10, 9};
//+
Line(46) = {9, 6};
//+
Line(47) = {11, 49};
//+
Line(48) = {49, 50};
//+
Line(49) = {50, 43};
//+
Line(50) = {43, 44};
//+
Line(51) = {44, 49};
//+
Line(52) = {49, 48};
//+
Line(53) = {48, 12};
//+
Line(54) = {48, 47};
//+
Line(55) = {47, 46};
//+
Line(56) = {46, 45};
//+
Line(57) = {45, 44};
//+
Line(58) = {48, 45};
//+
Line(59) = {45, 52};
//+
Line(60) = {52, 51};
//+
Line(61) = {51, 44};
//+
Line(62) = {10, 34};
//+
Line(63) = {34, 33};
//+
Line(64) = {33, 40};
//+
Line(65) = {40, 39};
//+
Line(66) = {39, 34};
//+
Line(67) = {34, 35};
//+
Line(68) = {35, 36};
//+
Line(69) = {36, 37};
//+
Line(70) = {37, 38};
//+
Line(71) = {38, 39};
//+
Line(72) = {39, 42};
//+
Line(73) = {42, 41};
//+
Line(74) = {41, 38};
//+
Line(75) = {38, 35};
//+
Line(76) = {35, 9};
//+
Line(77) = {2, 57};
//+
Line(78) = {57, 3};
//+
Line(79) = {3, 54};
//+
Line(80) = {54, 6};
//+
Line(81) = {6, 55};
//+
Line(82) = {55, 7};
//+
Line(83) = {7, 56};
//+
Line(84) = {56, 2};
//+
Line(85) = {57, 53};
//+
Line(86) = {53, 54};
//+
Line(87) = {53, 55};
//+
Line(88) = {53, 56};

//---------------------------------------------------------------------------
// This section describes the "Plane Surfaces", i.e., the 2D surfaces for meshing

Curve Loop(101) = {77, 85, 88, 84};
Plane Surface(101) = {101};

Curve Loop(102) = {78, 79, -86, -85};
Plane Surface(102) = {102};

Curve Loop(103) = {86, 80, 81, -87};
Plane Surface(103) = {103};

Curve Loop(104) = {-88, 87, 82, 83};
Plane Surface(104) = {104};

Curve Loop(2) = {42, 43, -78, -77, 41};
Plane Surface(2) = {2};

Curve Loop(3) = {6, 25, -15, -80, -79};
Plane Surface(3) = {3};

Curve Loop(4) = {-82, -81, -46, -45, 44};
Plane Surface(4) = {4};

Curve Loop(5) = {4, -84, -83, -13, -22};
Plane Surface(5) = {5};

Curve Loop(6) = {8, 27, -17, -26};
Plane Surface(6) = {6};

Curve Loop(7) = {39, 40, -8, 38};
Plane Surface(7) = {7};

Curve Loop(8) = {9, 28, -18, -27};
Plane Surface(8) = {8};

Curve Loop(9) = {17, -36, -35, 37};
Plane Surface(9) = {9};

Curve Loop(10) = {7, 26, -16, -25};
Plane Surface(10) = {10};

Curve Loop(11) = {67, -75, 71, 66};
Plane Surface(11) = {11};

Curve Loop(12) = {45, -76, -67, -62};
Plane Surface(12) = {12};

Curve Loop(13) = {68, 69, 70, 75};
Plane Surface(13) = {13};

Curve Loop(14) = {-71, -74, -73, -72};
Plane Surface(14) = {14};

Curve Loop(15) = {-63, -66, -65, -64};
Plane Surface(15) = {15};

Curve Loop(16) = {2, 21, -11, -20};
Plane Surface(16) = {16};

Curve Loop(17) = {33, 34, -2, 32};
Plane Surface(17) = {17};

Curve Loop(18) = {3, 22, -12, -21};
Plane Surface(18) = {18};

Curve Loop(19) = {11, -30, -29, 31};
Plane Surface(19) = {19};

Curve Loop(20) = {1, 20, -10, -19};
Plane Surface(20) = {20};

Curve Loop(21) = {-57, -58, -52, -51};
Plane Surface(21) = {21};

Curve Loop(22) = {-60, -59, 57, -61};
Plane Surface(22) = {22};

Curve Loop(23) = {-56, -55, -54, 58};
Plane Surface(23) = {23};

Curve Loop(24) = {52, 53, -42, 47};
Plane Surface(24) = {24};

Curve Loop(25) = {50, 51, 48, 49};
Plane Surface(25) = {25};

//---------------------------------------------------------------------------
// Creates physical surfaces. The number in () represents the material id.
// The expression list on the right hand side is the list
// of elementary surfaces created above.  This is what makes our mesh 2D.
Physical Surface(1) = {101, 102, 103, 104, 6, 11, 16, 21};
Physical Surface(2) = {2, 7, 12, 17, 22};
Physical Surface(3) = {3, 8, 13, 18, 23};
Physical Surface(4) = {4, 9, 14, 19, 24};
Physical Surface(5) = {5,10, 15, 20, 25};

//---------------------------------------------------------------------------
// Parameters for the meshing
Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 2.0;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
