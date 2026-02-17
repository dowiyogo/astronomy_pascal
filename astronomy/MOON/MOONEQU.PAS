(*-----------------------------------------------------------------------*)
(* MOONEQU: geocentric equatorial coordinates of the Moon                *)
(*          referred to the true equinox of date                         *)
(*   T   time in Julian centuries ephemeris time since J2000             *)
(*       ( T = (JD-2451545.0)/36525 )                                    *)
(*   RA  right ascension (deg)                                           *)
(*   DEC declination (deg)                                               *)
(*   R   distance (in earth radii)                                       *)
(*-----------------------------------------------------------------------*)
PROCEDURE MOONEQU(T:REAL;VAR RA,DEC,R:REAL);
  VAR L,B,X,Y,Z: REAL;
  BEGIN
    MOON(T,L,B,R);              (* ecliptic coordinates (mean equinox    *)
    CART(R,B,L,X,Y,Z);          (*                              of date) *)
    ECLEQU(T,X,Y,Z);            (* transform into equatorial coordinates *)
    NUTEQU(T,X,Y,Z);            (* nutation                              *)
    POLAR(X,Y,Z,R,DEC,RA);
  END;

