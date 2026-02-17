(*-----------------------------------------------------------------------*)
(* Unit MOOLIB: lunar orbit                                              *)
(*-----------------------------------------------------------------------*)

UNIT MOOLIB;

INTERFACE

  USES MATLIB, PNULIB, SPHLIB;

  PROCEDURE MINI_MOON  ( T: REAL; VAR RA,DEC: REAL );
  PROCEDURE MOON       ( T: REAL; VAR LAMBDA,BETA,R: REAL );
  PROCEDURE MOONEQU    ( T: REAL; VAR RA,DEC,R: REAL );
  PROCEDURE T_FIT_MOON ( TA,TB: REAL; N: INTEGER;
                         VAR RA_POLY,DE_POLY,R_POLY: TPOLYNOM );

IMPLEMENTATION

  (*$I \MOON\MINIMOON.PAS    *)
  (*$I \MOON\MOON.PAS        *)
  (*$I \MOON\MOONEQU.PAS     *)
  (*$I \DIVERSE\TFITMOON.PAS *)

END.
