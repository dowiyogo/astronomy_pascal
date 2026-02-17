(*-----------------------------------------------------------------------*)
(* T_FIT_SUN: approximates the equatorial coordinates of the             *)
(*            Sun by Chebyshev expansions for a given period of time     *)
(*                                                                       *)
(*  TA     : first date (in Julian centuries since J2000)                *)
(*  TB     : last date ( TB < TA + 1 year )                              *)
(*  N      : highest order                                               *)
(*  RA_POLY: coefficients for right ascension                            *)
(*  DE_POLY: coefficients for declination                                *)
(*  R_POLY : coefficients for geocentric distance                        *)
(*-----------------------------------------------------------------------*)
PROCEDURE T_FIT_SUN ( TA,TB: REAL; N: INTEGER;
                      VAR RA_POLY,DE_POLY,R_POLY: TPOLYNOM);
  BEGIN
    T_FIT_LBR (SUNEQU,TA,TB,N,RA_POLY,DE_POLY,R_POLY);
  END;
