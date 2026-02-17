(*-----------------------------------------------------------------------*)
(* Unit P69LIB: Saturn, Uranus, Neptune, Pluto                           *)
(*-----------------------------------------------------------------------*)

UNIT P69LIB;

INTERFACE

  PROCEDURE SAT200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE URA200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE NEP200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE PLU200 ( T: REAL; VAR L,B,R: REAL );

IMPLEMENTATION

  (*$I \PLANETS\SAT200.PAS *)
  (*$I \PLANETS\URA200.PAS *)
  (*$I \PLANETS\NEP200.PAS *)
  (*$I \PLANETS\PLU200.PAS *)

END.
