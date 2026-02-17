(*-----------------------------------------------------------------------*)
(* Unit PNULIB: precession and nutation                                  *)
(*-----------------------------------------------------------------------*)

UNIT PNULIB;

INTERFACE

  USES MATLIB;

  PROCEDURE NUTEQU    ( T: REAL; VAR X,Y,Z: REAL );
  PROCEDURE PMATECL   ( T1,T2: REAL; VAR A: REAL33 );
  PROCEDURE PMATEQU   ( T1,T2: REAL; VAR A: REAL33 );
  PROCEDURE PN_MATRIX ( T0,T: REAL; VAR A: REAL33 );
  PROCEDURE PRECART   ( A: REAL33; VAR X,Y,Z: REAL );

IMPLEMENTATION

  (*$I \PRECNUT\NUTEQU.PAS   *)
  (*$I \PRECNUT\PMATECL.PAS  *)
  (*$I \PRECNUT\PMATEQU.PAS  *)
  (*$I \PRECNUT\PNMATRIX.PAS *)
  (*$I \PRECNUT\PRECART.PAS  *)

END.
