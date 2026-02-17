(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Unit DELIB:                                                           *)
(*                                                                       *)
(*  Pascal implementation of the variable order variable stepsize        *)
(*  Adams-Bashforth-Moulton code DE for the numerical integration of     *)
(*  first order differential equations by Shampine and Gordon.           *)
(*                                                                       *)
(*  The original Fortran code is described in                            *)
(*                                                                       *)
(*  - Shampine L.F., Gordon M.K.: Computer solution of ordinary          *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*  - Shampine L.F., Gordon M.K.: Computer Loesung gewoehnlicher         *)
(*    Differentialgleichungen; Vieweg Verlag, Braunschweig (1984).       *)
(*  - Shampine L.F., Watts H.A.: DEPAC - Design of a user oriented       *)
(*    package of ODE solvers; SAND79-2374, Sandia Laboratories (1979).   *)
(*                                                                       *)
(*  The following functions and procedures are required to run DE:       *)
(*                                                                       *)
(*    INTRP           Interpolation routine                              *)
(*    STEP            Basic integration routine                          *)
(*    DE              User called driver routine                         *)
(*                                                                       *)
(*  Furthermore, various type and constant definitions are required that *)
(*  have to be included into the declaration part of the main program.   *)
(*                                                                       *)
(*  The constant DE_EQN specifies the maximum dimension of the system of *)
(*  ordinary differential equations that can be integrated. TWOU and     *)
(*  FOURU have to be defined as two and four times the accuracy of the   *)
(*  REAL data type.                                                      *)
(*                                                                       *)
(*    CONST DE_EQN = 6;  TWOU = 4.0E-12;  FOURU = 8.0E-12;               *)
(*                                                                       *)
(*  The DE_EQN_VECTOR type defines a vector of dimension DE_EQN for      *) 
(*  state variables, while DE_12_VECTOR, DE_13_VECTOR and DE_PHI_VECTOR  *)
(*  are array types used to define various auxiliary variables.          *)
(*  DE_WORKSPACE_RECORD finally defines a record structure used to hold  *)
(*  all data that have to be stored globally between successive calls of *)
(*  the DE procedure.                                                    *)
(*                                                                       *)
(*    TYPE  DE_EQN_VECTOR = ARRAY[1..DE_EQN] OF REAL;                    *)
(*          DE_12_VECTOR  = ARRAY[1..12] OF REAL;                        *)
(*          DE_13_VECTOR  = ARRAY[1..13] OF REAL;                        *)
(*          DE_PHI_ARRAY  = ARRAY[1..DE_EQN,1..16] OF REAL;              *)
(*          DE_WORKSPACE_RECORD                                          *)
(*                        = RECORD                                       *)
(*                            YY,WT,P,YP,YPOUT     :  DE_EQN_VECTOR;     *)
(*                            PHI                  :  DE_PHI_ARRAY;      *)
(*                            ALPHA,BETA,V,W,PSI   :  DE_12_VECTOR;      *)
(*                            SIG,G                :  DE_13_VECTOR;      *)
(*                            X,H,HOLD,TOLD,DELSGN :  REAL;              *)
(*                            NS,K,KOLD,ISNOLD     :  INTEGER;           *)
(*                            PHASE1,START,NORND   :  BOOLEAN;           *)
(*                          END;                                         *)
(*                                                                       *)
(*                                                                       *)
(*  The function to be integrated must be declared as a procedure of     *)
(*  the form PROCEDURE ( VAR X: REAL; VAR Y,YP: DE_EQN_VECTOR ); and     *)
(*  the same definition must be used in the procedure header of STEP     *)
(*  and DE.                                                              *)
(*                                                                       *)
(*  Using function EPSMACH the machine accuracy U=EPSMACH can be         *)
(*  computed to determine the proper values of constants TWOU and FOURU  *)
(*  for other compilers or floating point data types.                    *)
(*                                                                       *)
(*  Turbo Pascal adaption:                                               *)
(*                                                                       *)
(*  Turbo Pascal uses special data types to pass a function or procedure *)
(*  as formal argument of another function or procedure. For this        *)
(*  purpose the following procedure type is declared and used in the     *)
(*  specification of STEP and DE:                                        *)
(*                                                                       *)
(*      DE_FUNC = PROCEDURE ( X: REAL; VAR Y,YP: DE_EQN_VECTOR );        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)


UNIT DELIB;

INTERFACE

  (*---------------------------------------------------------------------*)
  (* Sample dimensions and type declarations for use with DE,STEP,INTRP  *)
  (*---------------------------------------------------------------------*)

  CONST DE_EQN = 6;           (* Maximum number of ODEs to be integrated *)
        TWOU   = 4.0E-12;     (* Turbo Pascal REAL data type accuracy    *)
        FOURU  = 8.0E-12;     (* u=2**(-39)=1.81E-12                     *)

  TYPE  DE_EQN_VECTOR       = ARRAY[1..DE_EQN] OF REAL;
        DE_FUNC = PROCEDURE ( X: REAL; Y: DE_EQN_VECTOR; VAR YP: DE_EQN_VECTOR );

        DE_12_VECTOR        = ARRAY[1..12] OF REAL;
        DE_13_VECTOR        = ARRAY[1..13] OF REAL;
        DE_PHI_ARRAY        = ARRAY[1..DE_EQN,1..16] OF REAL;
  
        DE_WORKSPACE_RECORD = RECORD 
                                YY,WT,P,YP,YPOUT     :  DE_EQN_VECTOR;
                                PHI                  :  DE_PHI_ARRAY;
                                ALPHA,BETA,V,W,PSI   :  DE_12_VECTOR;
                                SIG,G                :  DE_13_VECTOR;
                                X,H,HOLD,TOLD,DELSGN :  REAL;
                                NS,K,KOLD,ISNOLD     :  INTEGER;
                                PHASE1,START,NORND   :  BOOLEAN;
                              END;
  

  (*---------------------------------------------------------------------*)
  (* Procedures and Functions                                            *)
  (*---------------------------------------------------------------------*)

  FUNCTION EPSMACH: REAL;

  PROCEDURE DE ( F                  :  DE_FUNC;
                 NEQN               :  INTEGER;
                 VAR Y              :  DE_EQN_VECTOR;
                 VAR T,TOUT         :  REAL;
                 VAR RELERR,ABSERR  :  REAL;
                 VAR IFLAG          :  INTEGER;
                 VAR WORK           :  DE_WORKSPACE_RECORD        );
  


IMPLEMENTATION

  (*$I \DE\EPSMACH.PAS *)
  (*$I \DE\INTRP.PAS   *)
  (*$I \DE\STEP.PAS    *)
  (*$I \DE\DE.PAS      *)

END.
