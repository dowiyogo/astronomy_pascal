(*-----------------------------------------------------------------------*)
(* POLAR: conversion of cartesian coordinates (x,y,z)                    *)
(*        into polar coordinates (r,theta,phi)                           *)
(*        (theta in [-90 deg,+90 deg]; phi in [0 deg,+360 deg])          *)
(*-----------------------------------------------------------------------*)
PROCEDURE POLAR(X,Y,Z:REAL;VAR R,THETA,PHI:REAL);
  VAR RHO: REAL;
  BEGIN
    RHO:=X*X+Y*Y;  R:=SQRT(RHO+Z*Z);  
    PHI:=ATN2(Y,X); IF PHI<0 THEN PHI:=PHI+360.0;
    RHO:=SQRT(RHO); THETA:=ATN2(Z,RHO);
  END;
