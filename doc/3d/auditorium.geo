Include "3DMacro.geo";

SetFactory("OpenCASCADE");

angle = Pi/2;
radiusOut = 30;
crop = 18;
radiusIn = crop/100*radiusOut;
height = 10;

// Build geometry

extWall = newv;
Cylinder(extWall) = {0, 0, 0, 0, 0, height, radiusOut, angle};
inWall = newv;
Box(inWall) = {-radiusIn, -radiusIn, 0, 2*radiusIn, 2*radiusIn, height};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{inWall}; }
vWall = newv;
BooleanDifference(vWall) = { Volume{extWall}; Delete; }{ Volume{inWall}; Delete; };

// Main stairs

stairOffset= 15;
height = 5;
nStep = 10;
stepH = height/nStep;
stepL = 1;
i=0;

For i In {0:nStep-1}
	vStair = newv;
	Cylinder(vStair) = {0, 0, i*stepH, 0, 0, stepH, radiusOut, angle};
	vStairIn = newv;
	Cylinder(vStairIn) = {0, 0, i*stepH, 0, 0, stepH, stairOffset+i*stepL, angle};
	vStairReverse = newv;
	BooleanDifference(vStairReverse) = { Volume{vStair}; Delete; }{ Volume{vStairIn}; Delete; };
	vOld = vWall;
	vWall = newv;
	BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStairReverse}; Delete; };
EndFor

// Stage

stageH = 0.2;
stageW = 8;
stageL = 4;

vStage = newv;
Box(vStage) = {radiusIn, 0-stageW/2, 0, stageL, stageW, stageH};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{vStage}; }
vOld = vWall; vWall = newv;
BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStage}; Delete; };

// Desk

deskH = 1.1;
deskW = 6;
deskL = 1.5;

vDesk = newv;
Box(vDesk) = {radiusIn+0.5*stageL, -(deskW)/2, stageH, deskL, deskW, deskH};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{vDesk}; }
vOld = vWall; vWall = newv;
BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vDesk}; Delete; };

// Column

numCol = 4;
angleCol = angle/(numCol+0.2);
For i In {0:numCol-1}
	Col = newv;
	Cylinder(Col) = {radiusOut-2.5, 0, height, 0, 0, height, 1, 2*Pi};
	Rotate {{0, 0, 1}, {0, 0, 0}, -(i-numCol/2.+0.5)*angleCol+Pi/4} { Volume{Col};}
	vOld = vWall; vWall = newv;
	BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{Col}; Delete; };
EndFor

// Table

stageH = 5;
stageW = stageW;
stageL = 0.2;

vStage = newv;
Box(vStage) = {radiusIn, -stageW/2, 3, stageL, stageW, stageH};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{vStage}; }
vOld = vWall; vWall = newv;
BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStage}; Delete; };


// Small stairs

stairAngle = Pi/20;

For i In {0:nStep-1}
	vStair = newv;
	Cylinder(vStair) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL, stairAngle*stairOffset/(stairOffset+i*stepL)};
	vStairIn = newv;
	Cylinder(vStairIn) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL-stepL/3, angle};
	vStairReverse = newv;
	BooleanDifference(vStairReverse) = { Volume{vStair}; Delete; }{ Volume{vStairIn}; Delete; };
	Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4-0.5*stairAngle*stairOffset/(stairOffset+i*stepL)} { Volume{vStairReverse};}
	vOld = vWall;
	vWall = newv;
	BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStairReverse}; Delete; };

	vStair = newv;
	Cylinder(vStair) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL, stairAngle*stairOffset/(stairOffset+i*stepL)};
	vStairIn = newv;
	Cylinder(vStairIn) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL-stepL/3, angle};
	vStairReverse = newv;
	BooleanDifference(vStairReverse) = { Volume{vStair}; Delete; }{ Volume{vStairIn}; Delete; };
	Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2-stairAngle*stairOffset/(stairOffset+i*stepL)} { Volume{vStairReverse};}
	vOld = vWall;
	vWall = newv;
	BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStairReverse}; Delete; };

	vStair = newv;
	Cylinder(vStair) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL, stairAngle*stairOffset/(stairOffset+i*stepL)};
	vStairIn = newv;
	Cylinder(vStairIn) = {0, 0, i*stepH, 0, 0, stepH/2, stairOffset+i*stepL-stepL/3, angle};
	vStairReverse = newv;
	BooleanDifference(vStairReverse) = { Volume{vStair}; Delete; }{ Volume{vStairIn}; Delete; };
	Rotate {{0, 0, 1}, {0, 0, 0}, 0} { Volume{vStairReverse};}
	vOld = vWall;
	vWall = newv;
	BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vStairReverse}; Delete; };
EndFor
//+

For j In {0:nStep-2}
	For i In {1:(5+(j%1))}
		vBox = newv;
		Box(vBox) = {stairOffset+j*stepL+stepL/3, -stageL/20, (j+1)*stepH, stageW/25, stageW/10, stepH};
		Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3.6+(i+((j+1)%1))*(Pi/(2*15))-Pi/4+(j%1)*Pi/(4*15)} { Volume{vBox}; }
		vOld = vWall;
		vWall = newv;
		BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vBox}; Delete; };
	EndFor
EndFor

For j In {0:nStep-2}
	For i In {1:(5+(j%1))}
		vBox = newv;
		Box(vBox) = {stairOffset+j*stepL+stepL/3, -stageL/20, (j+1)*stepH, stageW/25, stageW/10, stepH};
		Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3.85+(i+((j+1)%1))*(Pi/(2*15))+(j%1)*Pi/(4*15)} { Volume{vBox}; }
		vOld = vWall;
		vWall = newv;
		BooleanDifference(vWall) = { Volume{vOld}; Delete; }{ Volume{vBox}; Delete; };
	EndFor
EndFor



