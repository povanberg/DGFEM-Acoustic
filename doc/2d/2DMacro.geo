Macro "HoleCircle"

  pCir[0] = newp;  Point(pCir[0]) = {x,y,0,LcHole};
  pCir[1] = newp;  Point(pCir[1]) = {x+r,y,0,LcHole};
  pCir[2] = newp;  Point(pCir[2]) = {x,y+r,0,LcHole};
  pCir[3] = newp;  Point(pCir[3]) = {x-r,y,0,LcHole};
  pCir[4] = newp;  Point(pCir[4]) = {x,y-r,0,LcHole};
  pCir[5] = pCir[1];

  For j In{0:3}
    cCir[j] = newl;  Circle(cCir[j]) = {pCir[j+1],pCir[0],pCir[j+2]};
  EndFor

  cHoles[] += {cCir[]};
  clCir = newll; Curve Loop(clCir) = {cCir[]};
  clHoles[] += {clCir};

Return

Macro "HoleLong"

  pLon[0] = newp;  Point(pLon[0]) = {x+w/2.0,y,0,LcHole};
  pLon[1] = newp;  Point(pLon[1]) = {x+w/2.0,y-r,0,LcHole};
  pLon[2] = newp;  Point(pLon[2]) = {x+w/2.0+r,y,0,LcHole};
  pLon[3] = newp;  Point(pLon[3]) = {x+w/2.0,y+r,0,LcHole};

  pLon[4] = newp;  Point(pLon[4]) = {x-w/2.0,y,0,LcHole};
  pLon[5] = newp;  Point(pLon[5]) = {x-w/2.0,y+r,0,LcHole};
  pLon[6] = newp;  Point(pLon[6]) = {x-w/2.0-r,y,0,LcHole};
  pLon[7] = newp;  Point(pLon[7]) = {x-w/2.0,y-r,0,LcHole};

  For j In{0:1}
    cLon[j] = newl;  Circle(cLon[j]) = {pLon[j+1],pLon[0],pLon[j+2]};
    cLon[j*3+2] = newl;  Line(cLon[j*3+2]) = {pLon[j*4+3],pLon[5-j*4]};
    cLon[j+3] = newl;  Circle(cLon[j+3]) = {pLon[j+5],pLon[4],pLon[j+6]};
  EndFor

  cHoles[] += {cLon[]};
  clLon = newll; Curve Loop(clLon) = {cLon[]};
  clHoles[] += {clLon};

Return
