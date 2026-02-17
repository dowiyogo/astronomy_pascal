(*-----------------------------------------------------------------------*)
(* ETMINUT: Difference of Ephemeris Time and Universal Time              *)
(*          (polynomial approximation valid from 1900 to 1995)           *)
(*          T:     time in Julian centuries since J2000                  *)
(*                 ( = (JD-2451545.0)/36525.0 )                          *)
(*          DTSEC: DT=ET-UT in sec (only if VALID=TRUE)                  *)
(*          VALID: TRUE for times between 1900 and 1995, FALSE otherwise *)
(*-----------------------------------------------------------------------*)
PROCEDURE ETMINUT(T: REAL; VAR DTSEC: REAL; VAR VALID: BOOLEAN);
  BEGIN
    VALID := ( (-1.0<=T) AND (T<=-0.05) );
    IF (VALID) THEN
      BEGIN
        DTSEC := ((((-449.50*T-783.42)*T-387.70)*T)+13.34)*T+62.14;
      END;
  END;

