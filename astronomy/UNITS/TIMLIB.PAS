(*-----------------------------------------------------------------------*)
(* Unit TIMLIB: time and calendar calculations                           *)
(*-----------------------------------------------------------------------*)
UNIT TIMLIB;

INTERFACE

  PROCEDURE CALDAT  ( MJD: REAL; VAR DAY,MONTH,YEAR: INTEGER;
                      VAR HOUR: REAL );
  PROCEDURE ETMINUT ( T:  REAL; VAR DTSEC:  REAL; VAR VALID:  BOOLEAN );
  FUNCTION LMST     ( MJD,LAMBDA: REAL ): REAL;
  FUNCTION MJD      ( DAY,MONTH,YEAR: INTEGER; HOUR: REAL ): REAL;

IMPLEMENTATION

  FUNCTION LONG_TRUNC ( X: REAL ): REAL;
    BEGIN  LONG_TRUNC := INT(X);  END;

  (*$I \TIME\CALDAT.PAS  *)
  (*$I \TIME\ETMINUT.PAS *)
  (*$I \TIME\LMST.PAS    *)
  (*$I \TIME\MJD.PAS     *)

END.
