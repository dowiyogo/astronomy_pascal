(*-----------------------------------------------------------------------*)
(* Unit P15LIB: Mercury, Venus, Mars, Jupiter                            *)
(*-----------------------------------------------------------------------*)

UNIT P15LIB;

INTERFACE

  USES MATLIB;

  PROCEDURE GEOCEN ( T, LP,BP,RP, LS,BS,RS: REAL; IPLAN,IMODE: INTEGER;
                     VAR XP,YP,ZP, XS,YS,ZS, X,Y,Z,DELTA0: REAL );
  PROCEDURE MER200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE VEN200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE MAR200 ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE JUP200 ( T: REAL; VAR L,B,R: REAL );

IMPLEMENTATION

  (*$I \PLANETS\GEOCEN.PAS *)
  (*$I \PLANETS\MER200.PAS *)
  (*$I \PLANETS\VEN200.PAS *)
  (*$I \PLANETS\MAR200.PAS *)
  (*$I \PLANETS\JUP200.PAS *)

END.
