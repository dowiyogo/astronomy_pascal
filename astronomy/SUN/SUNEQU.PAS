(*-----------------------------------------------------------------------*)
(* SUNEQU: apparent equatorial coordinates of the Sun                    *)
(*         (right ascension RA, declination DEC in deg, R in AU)         *)
(*         T in Julian centuries since J2000 (T=(JD-2451545.0)/36525)    *)
(*-----------------------------------------------------------------------*)
PROCEDURE SUNEQU(T:REAL;VAR RA,DEC,R:REAL);
  VAR DT,L,B,X,Y,Z: REAL;
  BEGIN
    DT := (8.32/1440.0)/36525.0; (* light-time correction of 8.32 min    *)
    SUN200(T-DT,L,B,R);          (* geocentric ecliptic coordinates      *)
    CART(R,B,L,X,Y,Z);           (* cartesian ecliptic coordinates       *)
    ECLEQU(T,X,Y,Z);             (* equatorial ecliptic coordinates      *)
    NUTEQU(T,X,Y,Z);             (* correction for nutation              *)
    POLAR(X,Y,Z,R,DEC,RA);       (* spherical coordinates;               *)
  END;                           (* true equinox of date                 *)

