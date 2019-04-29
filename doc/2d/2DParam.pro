DefineConstant[BcBord = {0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "Parameters/1BC"}];
DefineConstant[LcBord = {1.0, Min 1e-2, Max 10, Highlight "Green", Name "Parameters/2Mesh size", Units "m"}];
DefineConstant[xBord = {30.0, Min 0, Max 100, Highlight "Green", Name "Parameters/3Width", Units "m"}];
DefineConstant[yBord = {15.0, Min 0, Max 100, Highlight "Green", Name "Parameters/4Height", Units "m"}];

DefineConstant[BcHole = {0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "Parameters/Holes/1BC"}];
DefineConstant[LcHole = {1.0, Min 1e-2, Max 100, Highlight "Green", Name "Parameters/Holes/2Mesh size", Units "m"}];
DefineConstant[nHoleCircle = {0, Min 0, Max 20, Step 1, Highlight "Green", Name "Parameters/Holes/3Circle", Units "-"}];
DefineConstant[nHoleLong = {0, Min 0, Max 20, Step 1, Highlight "Green", Name "Parameters/Holes/3Long", Units "-"}];

For i In {1:nHoleCircle}
  DefineConstant[
    XC~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Circle %g/1x", i), Units "m"},
    YC~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Circle %g/2y", i), Units "m"},
    ZC~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Circle %g/3z", i), Units "m"},
    RC~{i} = {1.0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Circle %g/4r", i), Units "m"}
  ];
EndFor

For i In {1:nHoleLong}
  DefineConstant[
    XL~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Long %g/1x", i), Units "m"},
    YL~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Long %g/2y", i), Units "m"},
    ZL~{i} = {0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Long %g/3z", i), Units "m"},
    RL~{i} = {1.0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Long %g/4r", i), Units "m"},
    WL~{i} = {1.0, Min 0, Max 100, Highlight "Green", Name Sprintf("Parameters/Holes/Long %g/4w", i), Units "m"}
  ];
EndFor
