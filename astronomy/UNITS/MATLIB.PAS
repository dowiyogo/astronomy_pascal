(*-----------------------------------------------------------------------*)
(* Unit MATLIB : mathematical functions and subroutines                  *)
(*-----------------------------------------------------------------------*)

UNIT MATLIB;

INTERFACE

  (* vectors and matrices *)

  TYPE INDEX  = (X,Y,Z);
       VECTOR = ARRAY[INDEX] OF REAL;
       REAL3  = ARRAY[1.. 3] OF REAL;
       REAL33 = ARRAY[1.. 3] OF REAL3;
       MAT3X  = ARRAY[1.. 3] OF VECTOR;

  (* Chebyshev polynomials *)

  CONST MAX_TP_DEG  = 13;                        (* maximum degree       *)
  TYPE  TPOLYNOM = RECORD                        (* Chebyshev polynomial *)
                     M  : INTEGER;                       (* degree       *)
                     A,B: REAL;                          (* interval     *)
                     C  : ARRAY [0..MAX_TP_DEG] OF REAL; (* coefficients *)
                   END;

  (* Vector and matrix for least squares systems:                    *)
  (* minimum dimensions for use with program FOTO are given below;   *)
  (* they may be increased as needed by the user.                    *)

  TYPE  LSQVEC =  ARRAY[1..3] OF REAL;
        LSQMAT =  ARRAY[1..30,1..5] OF REAL;


  (* functions and procedures *)

  FUNCTION  ACS    ( X: REAL ): REAL;
  FUNCTION  ASN    ( X: REAL ): REAL;
  FUNCTION  ATN    ( X: REAL ): REAL;
  FUNCTION  ATN2   ( Y,X: REAL ): REAL;
  PROCEDURE CART   ( R,THETA,PHI: REAL; VAR X,Y,Z: REAL );
  PROCEDURE CROSS  ( A,B: VECTOR; VAR C: VECTOR );
  FUNCTION  CS     ( X: REAL ): REAL;
  FUNCTION  CUBR   ( X: REAL ): REAL;
  FUNCTION  DOT    ( A,B: VECTOR ): REAL;
  PROCEDURE DDD    ( D,M: INTEGER;S: REAL;VAR DD: REAL );
  PROCEDURE DMS    ( DD: REAL;VAR D,M: INTEGER;VAR S: REAL );
  PROCEDURE LSQFIT ( A: LSQMAT; N,M: INTEGER; VAR S: LSQVEC );
  FUNCTION  NORM   ( A: VECTOR ): REAL;
  PROCEDURE POLAR  ( X,Y,Z: REAL; VAR R,THETA,PHI: REAL );
  PROCEDURE QUAD   ( Y_MINUS,Y_0,Y_PLUS:  REAL;
                     VAR XE,YE,ZERO1,ZERO2: REAL; VAR NZ: INTEGER );
  FUNCTION  SN     ( X: REAL ): REAL;
  FUNCTION  T_EVAL ( F:  TPOLYNOM; X:  REAL ):  REAL;
  FUNCTION  TN     ( X: REAL ): REAL;

IMPLEMENTATION

  (*$I \MATH\ACS.PAS    *)
  (*$I \MATH\ASN.PAS    *)
  (*$I \MATH\ATN.PAS    *)
  (*$I \MATH\ATN2.PAS   *)
  (*$I \MATH\CART.PAS   *)
  (*$I \MATH\CROSS.PAS  *)
  (*$I \MATH\CS.PAS     *)
  (*$I \MATH\CUBR.PAS   *)
  (*$I \MATH\DOT.PAS    *)
  (*$I \MATH\DDD.PAS    *)
  (*$I \MATH\DMS.PAS    *)
  (*$I \MATH\LSQFIT.PAS *)
  (*$I \MATH\NORM.PAS   *)
  (*$I \MATH\POLAR.PAS  *)
  (*$I \MATH\QUAD.PAS   *)
  (*$I \MATH\SN.PAS     *)
  (*$I \MATH\TEVAL.PAS  *)
  (*$I \MATH\TN.PAS     *)

END.
