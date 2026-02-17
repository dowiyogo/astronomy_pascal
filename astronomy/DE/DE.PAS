(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* DE                                                                    *)
(*                                                                       *)
(* Driver routine for the variable order and step size Adams-Bashforth-  *)
(* Moulton method of Shampine and Gordon.                                *)
(*                                                                       *)
(*   F       Function to be integrated                                   *)
(*   NEQN    Number of differential equations                            *)
(*   Y       State vector                                                *)
(*   T       Time                                                        *)
(*   TOUT    Output time                                                 *)
(*   RELERR  Relative tolerance                                          *)
(*   ABSERR  Absolute tolerance                                          *)
(*   IFLAG   Return code                                                 *)
(*   WORK    Work area for global storage between subsequent calls       *)
(*                                                                       *)
(* Note:                                                                 *)
(*                                                                       *)
(* Subroutine DE is used for the numerical integration of an ordinary    *)
(* differential equation y'=f(x,y) of first order, which must be         *)
(* provided as a subroutine of the form                                  *)
(*   PROCEDURE F ( X: REAL; Y: DE_EQN_VECTOR; VAR YP: DE_EQN_VECTOR );   *)
(* Aside from DE itself the subroutines STEP and INTRP as well as some   *)
(* global constant and type definitions are required, which are          *)
(* described in detail in DELIB.PAS.                                     *)
(*                                                                       *)
(* On the first call of DE the initial time is passed in T, while Y      *)
(* contains the corresponding state vector. Furthermore the desired      *)
(* output time TOUT and the required relative and absolute accuracy      *)
(* RELERR and ABSERR need to be specified. The status variable IFLAG     *)
(* is set to 1 to inform DE that the integration of a new problem is     *)
(* started.                                                              *)
(*                                                                       *)
(* In normal cases the variable Y subsequently contains the state vector *)
(* at time TOUT. Simultaneously T is assigned the value TOUT, while the  *)
(* status flag IFLAg is set to 2 (successful step). In all continuation  *)
(* steps it is only required for the user to specify a new value TOUT,   *)
(* while all remaining variables (most notably IFLAG) are left           *)
(* unchanged. IFLAG needs only be reset to one in case of a new problem  *)
(* or when changing the direction of integration.                        *)
(*                                                                       *)
(* If DE returns a value of IFLAG different from 2, one of the following *)
(* cases may have occured:                                               *)
(*                                                                       *)
(*   IFLAG = 3  TOUT has not been reached, since the requested           *)
(*              tolerance RELERR or ABSERR was too small. Both           *)
(*              values have been increased for subsequent                *)
(*              computations.                                            *)
(*   IFLAG = 4  TOUT has not been reached, since more than               *)
(*              MAXNUM=500 steps have been required internally.          *)
(*   IFLAG = 5  TOUT has not been reached, since the differential        *)
(*              equation appears to be stiff. This should not happen     *)
(*              when applying DE to problems in celestial mechanics.     *)
(*   IFLAG = 6  Illegal input parameter (e.g. T=TOUT)                    *)
(*                                                                       *)
(* In all cases except IFLAG=6 the integration may simply be continued   *)
(* by calling DE again without changing any parmeter. The interruption   *)
(* is mainly intended to call the user's attention to possible problems  *)
(* that may have occured and e.g. avoid infinite integrations.           *)
(*                                                                       *)
(* The work array WORK is used to store various quantities between       *)
(* subsequent calls of DE and is not changed by the user.                *)
(*                                                                       *)
(* Note that "F:DE_FUNC;" in the specification of DE has to be replaced  *)
(* by "PROCEDURE F(X:REAL;Y:DE_EQN_VECTOR;VAR YP:DE_EQN_VECTOR);" for    *)
(* use with standard Pascal.                                             *)
(*                                                                       *)
(* A detailed description of the method and its Fortran implementation   *)
(* is given in                                                           *)
(*                                                                       *)
(*  - Shampine L.F., Gordon M.K.: Computer solution of ordinary          *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*  - Shampine L.F., Watts H.A.: DEPAC - Design of a user oriented       *)
(*    package of ODE solvers; SAND79-2374, Sandia Laboratories (1979).   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE DE ( F                  :  DE_FUNC;            (* Turbo Pascal *)
               NEQN               :  INTEGER;
               VAR Y              :  DE_EQN_VECTOR;
               VAR T,TOUT         :  REAL;
               VAR RELERR,ABSERR  :  REAL;
               VAR IFLAG          :  INTEGER;
               VAR WORK           :  DE_WORKSPACE_RECORD );

  LABEL 99;     (* Error exit at end of procedure *)

   
  CONST MAXNUM = 500;  (* Maximum number of steps *)


  VAR STIFF, CRASH                           :  BOOLEAN;
      ISN, NOSTEP, KLE4, L                   :  INTEGER;
      RELEPS, ABSEPS, TEND, ABSDEL, DEL, EPS :  REAL;


  (* SIGN: computes sign(B)*abs(A) *)

  FUNCTION SIGN(A,B: REAL): REAL; 
    BEGIN IF (B>=0.0) THEN SIGN:=ABS(A) ELSE SIGN:=-ABS(A); END;

  (* MIN: computes the minimum of A and B *)

  FUNCTION MIN(A,B: REAL): REAL;
    BEGIN IF A<B THEN MIN:=A ELSE MIN:=B; END;

  (* MAX: computes the maximum of A and B *)

  FUNCTION MAX(A,B: REAL): REAL;
    BEGIN IF A>B THEN MAX:=A ELSE MAX:=B; END;

  BEGIN

    WITH WORK DO   (* Use short form of record component name *)

      BEGIN
  
        (* Test for invalid parameters *)
        
        EPS   := MAX(RELERR,ABSERR);
        ISN   := TRUNC ( SIGN(1.1,IFLAG) );  (* should be +/- 1 *)
        IFLAG := ABS(IFLAG);
        
        IF ( ( NEQN   < 1   ) OR ( T      = TOUT ) OR
             ( RELERR < 0.0 ) OR ( ABSERR < 0.0  ) OR ( EPS = 0.0 ) OR
             ( IFLAG  = 0   ) OR ( IFLAG  > 5    ) OR
             ( (IFLAG<>1) AND (T<>TOLD) )               )  THEN 

          BEGIN  
        
            IFLAG:=6;  (* Set error code *)
            GOTO 99;   (* Exit           *)
        
          END;
        
        
        (* Set interval limits and step counter *)
        
        DEL    := TOUT - T; 
        ABSDEL := ABS(DEL);
        TEND   := T + 10.0*DEL;
        IF (ISN<0) THEN  TEND:=TOUT;
        NOSTEP := 0;
        KLE4   := 0;
        STIFF  := FALSE;
        RELEPS := RELERR/EPS;
        ABSEPS := ABSERR/EPS;
        
        IF ( (IFLAG=1) OR (ISNOLD<0) OR (DELSGN*DEL<=0.0) ) THEN  
        
          BEGIN 
        
            (* Set independent and dependent variables X and YY[*] for  *)
            (* steps. Set sign of integration direction. Initialize the *)
            (* step size.                                               *)
        
            START:=TRUE;
            X := T;
            FOR L:=1 TO NEQN DO YY[L]:=Y[L];
            DELSGN := SIGN(1.0,DEL);
            H := SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X);
        
          END;
        

        REPEAT 
        
          (* If already past output point, interpolate and return *)
          
          IF (ABS(X-T)>=ABSDEL) THEN
            BEGIN
              INTRP(X,YY,TOUT,Y,YPOUT,NEQN,KOLD,PHI,PSI);
              IFLAG := 2;
              T := TOUT;
              TOLD := T;
              ISNOLD := ISN;
              GOTO 99;     (* Exit *)
             END;
          
          (* If cannot go past TSTOP and sufficiently close, *)
          (* extrapolate and return                          *)
          
          IF ( (ISN<=0) AND (ABS(TOUT-X)<FOURU*ABS(X)) ) THEN
            BEGIN
              H := TOUT - X;
              F(X,YY,YP);
              FOR L:=1 TO NEQN DO Y[L]:=YY[L]+H*YP[L];
              IFLAG := 2;
              T := TOUT;
              TOLD := T;
              ISNOLD := ISN;
              GOTO 99;       (* Exit *)
            END;
          
          (* Monitor number of steps attempted *)
          
          IF (NOSTEP>=MAXNUM) THEN  
            (* a large amount of work has been expended *)
            BEGIN
              IFLAG := ISN*4;
              IF (STIFF) THEN  IFLAG:=ISN*5;
              FOR L:=1 TO NEQN DO  Y[L]:=YY[L];
              T := X;
              TOLD := T;
              ISNOLD := 1;
              GOTO 99;     (* Exit *)
            END;
          
          (* Limit step size, set weight vector and take a step *)
          
          H := SIGN(MIN(ABS(H),ABS(TEND-X)),H);
          FOR L:=1 TO NEQN DO  WT[L]:=RELEPS*ABS(YY[L])+ABSEPS;
          
          STEP ( X,YY,F,NEQN,H,EPS,WT,START,          
                 HOLD,K,KOLD,CRASH,PHI,P,YP,PSI,
                 ALPHA,BETA,SIG,V,W,G,PHASE1,NS,NORND );
          
          IF (CRASH) THEN   (* Tolerances too small *)
            BEGIN
              IFLAG   := ISN*3;
              RELERR  := EPS*RELEPS;
              ABSERR  := EPS*ABSEPS;
              FOR L:=1 TO NEQN DO Y[L]:=YY[L];
              T := X;
              TOLD   := T;
              ISNOLD := 1;
              GOTO 99;       (* Exit *)
            END;
          
          (* Stiffness test:                                  *)
          (* count number of consecutive steps taken with the *)
          (* order of the method being less or equal to four  *)
          
          NOSTEP := NOSTEP+1;
          KLE4   := KLE4+1;
          IF ( KOLD >   4 ) THEN KLE4:=0;
          IF ( KLE4 >= 50 ) THEN STIFF:=TRUE;
          
        UNTIL FALSE;
      

      END; (* with WORK do *) 

  99:

  END;

