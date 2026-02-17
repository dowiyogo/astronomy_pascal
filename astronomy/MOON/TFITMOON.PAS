(*-----------------------------------------------------------------------*)
(* T_FIT_MOON: approximates the equatorial coordinates                   *)
(*             of the Moon by Chebyshev expansions for a                 *)
(*             given period of time of at most one month                 *)
(*                                                                       *)
(*  TA     : first date (in Julian centuries since J2000)                *)
(*  TB     : last date ( TB < TA + 1 month )                             *)
(*  N      : highest order                                               *)
(*  RA_POLY: coefficients for right ascension                            *)
(*  DE_POLY: coefficients for declination                                *)
(*  R_POLY : coefficients for geocentric distance                        *)
(*-----------------------------------------------------------------------*)
PROCEDURE T_FIT_MOON ( TA,TB: REAL; N: INTEGER;
                       VAR RA_POLY,DE_POLY,R_POLY: TPOLYNOM);
  BEGIN
    T_FIT_LBR ( MOONEQU, TA,TB,N, RA_POLY,DE_POLY,R_POLY );
  END;

