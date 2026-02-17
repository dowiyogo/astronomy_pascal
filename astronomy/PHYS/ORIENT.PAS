(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* ORIENT:  returns the elements describing the planetocentric           *)
(*          coordinate system of a planet                                *)
(*                                                                       *)
(*   PLANET   Name of the planet                                         *)
(*   SYSTEM   System of rotation (I, II or III)                          *)
(*   T        Time in Julian centuries since J2000  (ET or TDB/TDT)      *)
(*   A        Right ascension of the axis of rotation                    *)
(*   D        Declination of the axis of rotation                        *)
(*   W        Orientation of the prime meridian with respect to the      *)
(*            intersection of the Earth's equator of J2000 and the       *)
(*            planetary equator of date                                  *)
(*   SENSE    Sense of rotation (DIRECT or RETROGRADE)                   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE ORIENT ( PLANET: PLANET_TYPE; SYSTEM: SYSTEM_TYPE; T: REAL;
                   VAR A,D,W: REAL; VAR SENSE: ROTATION_TYPE );

  VAR TD : REAL;

  BEGIN

    (* Compute right ascension and declination of the axis of    *)
    (* rotation with respect to the equator and equinox of J2000 *)

    CASE PLANET OF
      MERCURY : BEGIN A := 281.02 - 0.033*T; D :=  61.45 -0.005*T END;
      VENUS   : BEGIN A := 272.78          ; D :=  67.21          END;
      EARTH   : BEGIN A :=   0.00 - 0.641*T; D :=  90.00 -0.557*T END;
      MARS    : BEGIN A := 317.681- 0.108*T; D :=  52.886-0.061*T END;
      JUPITER : BEGIN A := 268.05 - 0.009*T; D :=  64.49 +0.003*T END;
      SATURN  : BEGIN A :=  40.66 - 0.036*T; D :=  83.52 -0.004*T END;
      URANUS  : BEGIN A := 257.43          ; D := -15.10          END;
      NEPTUNE : BEGIN A := 295.33          ; D :=  40.65          END;
      PLUTO   : BEGIN A := 311.63          ; D :=   4.18          END
    END;

    (* Compute orientation of the prime meridian *)

    TD := 36525.0 * T;
    CASE PLANET OF
      MERCURY : W := 329.710 +   6.1385025*TD;
      VENUS   : W := 159.910 -   1.4814205*TD;
      EARTH   : W := 100.21  + 360.9856123*TD;
      MARS    : W := 176.655 + 350.8919830*TD;
      JUPITER : CASE SYSTEM OF
                  SYS_I   : W :=  67.10 + 877.900*TD;
                  SYS_II  : W :=  43.30 + 870.270*TD;
                  SYS_III : W := 284.95 + 870.536*TD
                END;
      SATURN  : CASE SYSTEM OF
                  SYS_I, SYS_II : W := 227.2037 + 844.3000000*TD;
                  SYS_III       : W :=  38.90   + 810.7939024*TD
                END;
      URANUS  : W := 261.620 - 554.913*TD;  (* System III *)
      NEPTUNE : W := 107.210 + 468.75 *TD;  (* System III *)
      PLUTO   : W := 252.660 -  56.364*TD
    END;
    W := W/360.0;  W := 360.0*(W-TRUNC(W));

    (* Define sense of rotation *)

    CASE PLANET OF
      MERCURY,EARTH,MARS,
      JUPITER,SATURN,NEPTUNE :  SENSE := DIRECT;
      VENUS,URANUS,PLUTO     :  SENSE := RETROGRADE;
    END;

  END;

