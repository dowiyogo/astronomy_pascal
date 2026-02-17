(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* BRIGHT: Computes the apparent magnitude of a planet                   *)
(*                                                                       *)
(*  PLANET     Name of the planet                                        *)
(*  R          Heliocentric distance of the planet (AU)                  *)
(*  DELTA      Distance of the planet from the observer (AU)             *)
(*  PHI        Phase angle (deg)                                         *)
(*  DEC        Planetocentric latitude of the observer (deg)             *)
(*  DLONG      Difference of the planetocentric longitudes of the Sun    *)
(*             and the observer (deg)                                    *)
(*                                                                       *)
(* Magnitudes V(1,0) from Astronomical Almanac 1984. DEC and DLONG are   *)
(* only required to compute the apparent brightness of Saturn, which     *)
(* depends on the ring orientation.                                      *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

FUNCTION BRIGHT ( PLANET: PLANET_TYPE; R,DELTA,PHI,DEC,DLONG: REAL ): REAL;

  CONST LN10 = 2.302585093;    (* Natural logarithm of 10 *)

  VAR   P,SD,DL,MAG :  REAL;

  BEGIN

    P := PHI/100.0;

    CASE PLANET OF
      MERCURY : MAG := -0.42 + ( 3.80 - ( 2.73 - 2.00*P ) * P ) * P;
      VENUS   : MAG := -4.40 + ( 0.09 + ( 2.39 - 0.65*P ) * P ) * P;
      EARTH   : MAG := -3.86;
      MARS    : MAG := -1.52 + 1.6 * P;
      JUPITER : MAG := -9.40 + 0.5 * P;
      SATURN  : BEGIN
                  SD := ABS ( SN(DEC) );
                  DL := ABS(DLONG/100.0); IF (DL>1.8) THEN DL:=ABS(DL-3.6);
                  MAG := -8.88 - 2.60 * SD + 1.25 * SD*SD + 4.40 * DL;
                END;
      URANUS  : MAG := -7.19;
      NEPTUNE : MAG := -6.87;
      PLUTO   : MAG := -1.0;
    END;

    BRIGHT := MAG + 5.0 * LN ( R*DELTA ) / LN10;

  END;


