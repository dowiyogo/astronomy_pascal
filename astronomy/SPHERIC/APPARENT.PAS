(*-----------------------------------------------------------------------*)
(* APPARENT  : apparent coordinates of a star                            *)
(*   PNMAT   : matrix for precession and nutation                        *)
(*   VX,VY,VZ: velocity of the earth (equatorial coord.; in units of c)  *)
(*   RA,DEC  : right ascension and declination                           *)
(*-----------------------------------------------------------------------*)
PROCEDURE APPARENT ( PNMAT: REAL33; VX,VY,VZ: REAL; VAR RA,DEC: REAL );
  VAR X,Y,Z,R: REAL;
  BEGIN
    CART ( 1.0, DEC,RA, X,Y,Z ); (* cartesian coordinates of the star    *)
    PRECART ( PNMAT , X,Y,Z );   (* correct for precession and nutation  *)
    X:=X+VX; Y:=Y+VY; Z:= Z+VZ;  (* aberration                           *)
    POLAR( X,Y,Z, R,DEC,RA );    (* apparent right ascension,declination *)
  END;

