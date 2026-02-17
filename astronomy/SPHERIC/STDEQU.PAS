(*-----------------------------------------------------------------------*)
(* STDEQU: transformation from standard coordinates into                 *)
(*         equatorial coordinates                                        *)
(*   RA0,DEC0: right ascension and declination of the optical axis (deg) *)
(*   XX,YY:    standard coordinates                                      *)
(*   RA,DEC:   right ascension and declination (deg)                     *)
(*-----------------------------------------------------------------------*)
PROCEDURE STDEQU ( RA0,DEC0,XX,YY: REAL; VAR RA,DEC: REAL);
  BEGIN
    RA  := RA0 + ATN ( -XX / (CS(DEC0)-YY*SN(DEC0)) );
    DEC := ASN ( (SN(DEC0)+YY*CS(DEC0))/SQRT(1.0+XX*XX+YY*YY) );
  END;

