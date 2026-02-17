(*-----------------------------------------------------------------------*)
(* EQUSTD: transformation of equatorial coordinates into                 *)
(*         standard coordinates                                          *)
(*   RA0,DEC0: right ascension and declination of the optical axis (deg) *)
(*   RA,DEC:   right ascension and declination (deg)                     *)
(*   XX,YY:    standard coordinates                                      *)
(*-----------------------------------------------------------------------*)
PROCEDURE EQUSTD ( RA0,DEC0,RA,DEC: REAL; VAR XX,YY: REAL);
  VAR C: REAL;
  BEGIN
    C  := CS(DEC0)*CS(DEC)*CS(RA-RA0)+SN(DEC0)*SN(DEC);
    XX := - ( CS(DEC)*SN(RA-RA0) ) / C;
    YY := - ( SN(DEC0)*CS(DEC)*CS(RA-RA0)-CS(DEC0)*SN(DEC) ) / C;
  END;

