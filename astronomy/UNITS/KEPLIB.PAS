(*-----------------------------------------------------------------------*)
(* Unit KEPLIB: functions and subroutines for Keplerian orbits           *)
(*-----------------------------------------------------------------------*)

UNIT KEPLIB;

INTERFACE

  USES MATLIB, SPHLIB;

  FUNCTION  ECCANOM  ( MAN,ECC: REAL ): REAL;
  PROCEDURE ELEMENT  ( JDA,JDB: REAL; RA,RB: VECTOR;
                       VAR TP,Q,ECC,INC,LAN,AOP: REAL );
  PROCEDURE ELLIP    ( M,A,ECC: REAL; VAR X,Y,VX,VY: REAL );
  FUNCTION  FIND_ETA ( RA,RB: VECTOR; TAU: REAL ): REAL;
  FUNCTION  HYPANOM  ( MH,ECC: REAL ): REAL;
  PROCEDURE HYPERB   ( T0,T,A,ECC: REAL; VAR X,Y,VX,VY: REAL );
  PROCEDURE KEPLER   ( T0,T,Q,ECC: REAL; PQR: REAL33;
                       VAR X,Y,Z,VX,VY,VZ: REAL );
  PROCEDURE PARAB    ( T0,T,Q,ECC: REAL; VAR X,Y,VX,VY: REAL );
  PROCEDURE STUMPFF  ( E2: REAL; VAR C1,C2,C3: REAL );
  PROCEDURE XYZKEP   ( X,Y,Z, VX,VY,VZ           :  REAL;
                       VAR AX,ECC,INC,LAN,AOP,MA :  REAL );

IMPLEMENTATION

  (*$I \KEPLER\ECCANOM.PAS *)
  (*$I \KEPLER\ELEMENT.PAS *)
  (*$I \KEPLER\ELLIP.PAS   *)
  (*$I \KEPLER\FINDETA.PAS *)
  (*$I \KEPLER\HYPANOM.PAS *)
  (*$I \KEPLER\HYPERB.PAS  *)
  (*$I \KEPLER\KEPLER.PAS  *)
  (*$I \KEPLER\PARAB.PAS   *)
  (*$I \KEPLER\STUMPFF.PAS *)
  (*$I \KEPLER\XYZKEP.PAS  *)

END.
