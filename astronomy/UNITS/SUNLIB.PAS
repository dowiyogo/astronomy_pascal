(*-----------------------------------------------------------------------*)
(* Unit SUNLIB: solar orbit                                              *)
(*-----------------------------------------------------------------------*)

UNIT SUNLIB;

INTERFACE

  USES MATLIB, PNULIB, SPHLIB;

  PROCEDURE MINI_SUN  ( T: REAL; VAR RA,DEC: REAL );
  PROCEDURE SUN200    ( T: REAL; VAR L,B,R: REAL );
  PROCEDURE SUNEQU    ( T: REAL; VAR RA,DEC,R: REAL );
  PROCEDURE T_FIT_SUN ( TA,TB: REAL; N: INTEGER;
                        VAR RA_POLY,DE_POLY,R_POLY: TPOLYNOM );

IMPLEMENTATION

  (*$I \SUN\MINISUN.PAS     *)
  (*$I \SUN\SUN200.PAS      *)
  (*$I \SUN\SUNEQU.PAS      *)
  (*$I \DIVERSE\TFITSUN.PAS *)
     
END.
