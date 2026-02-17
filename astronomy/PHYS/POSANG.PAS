(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* POSANG: computes the position angle from given directions             *)
(*                                                                       *)
(*   X,Y,Z      Coordinates of the planet relative to the observer       *)
(*   DX,DY,DZ   Direction vector                                         *)
(*   POSANG     Position angle (0<=POSANG<360deg)                        *)
(*                                                                       *)
(* Both vectors must be given in a common coordinate system              *)
(* (e.g. mean equator and equinox of date)                               *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
FUNCTION POSANG ( X,Y,Z, DX,DY,DZ : REAL ): REAL;
  VAR C,S,PHI : REAL;
  BEGIN
    C := ( (-X*Z)*DX + (-Y*Z)*DY + (X*X+Y*Y)*DZ ) / SQRT (X*X+Y*Y+Z*Z);
    S := ( (-Y)*DX + (X)*DY );
    PHI := ATN2 ( S, C );
    IF (PHI<0.0) THEN POSANG:=PHI+360.0 ELSE POSANG:=PHI
  END; 

