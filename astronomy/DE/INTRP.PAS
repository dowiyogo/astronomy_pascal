(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* INTRP: Interpolation routine for use with DE and STEP                 *)
(*                                                                       *)
(*   X      Independent variable                                         *)
(*   Y      Solution at X                                                *)
(*   XOUT   Independent variable for which the solution is requested     *)
(*   YOUT   Interpolated solution at XOUT                                *)
(*   YPOUT  Derivative at XOUT                                           *)
(*   NEQN   Number of differential equations                             *)
(*   KOLD   Auxiliary quantity (most recently used order)                *)
(*   PHI    Auxiliary quantity                                           *)
(*   PSI    Auxiliary quantity                                           *)
(*                                                                       *)
(* STEP approximates the solution of the differential equation near X by *)
(* a polynomial, which is evaluated by INTRP.                            *)
(*                                                                       *)
(* A detailed description of the method is given in                      *)
(*                                                                       *)
(*    L. F. Shampine, M. K. Gordon; Computer Solution of ordinary        *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE INTRP ( X              :  REAL;
                  Y              :  DE_EQN_VECTOR;
                  XOUT           :  REAL;
                  VAR YOUT,YPOUT :  DE_EQN_VECTOR;
                  NEQN,KOLD      :  INTEGER;
                  PHI            :  DE_PHI_ARRAY;
                  PSI            :  DE_12_VECTOR  );

  VAR I,J,JM1,KI,KIP1,L,LIMIT1 :  INTEGER;
      ETA,GAMMA,HI,PSIJM1      :  REAL;
      TEMP1,TEMP2,TEMP3,TERM   :  REAL;
      G,RHO,W                  :  DE_13_VECTOR;


  BEGIN

    G[1]   := 1.0;
    RHO[1] := 1.0;

    HI := XOUT - X;
    KI := KOLD + 1;
    KIP1 := KI + 1;

    (* Initialize W[*] for computing G[*] *)

    FOR I:=1 TO KI DO 
      BEGIN 
        TEMP1 := I;
        W[I] := 1.0/TEMP1;
      END;
    TERM := 0.0;

    (* compute G[*] *)

    FOR J:=2 TO KI DO 
      BEGIN
        JM1 := J - 1;
        PSIJM1 := PSI[JM1];
        GAMMA := (HI + TERM )/PSIJM1;
        ETA := HI/PSIJM1;
        LIMIT1 := KIP1 - J;
        FOR I:=1 TO LIMIT1 DO  W[I] := GAMMA*W[I] - ETA*W[I+1];
        G[J] := W[1];
        RHO[J] := GAMMA*RHO[JM1];
        TERM := PSIJM1;
      END;

    (* Interpolate for the solution YOUT and for *)
    (* the derivative of the solution YPOUT      *)

    FOR L:=1 TO NEQN DO 
      BEGIN
        YPOUT[L] := 0.0;
        YOUT[L]  := 0.0;
      END;

    FOR J:=1 TO KI DO
      BEGIN
        I := KIP1 - J;
        TEMP2 := G[I];
        TEMP3 := RHO[I];
        FOR L:=1 TO NEQN DO
          BEGIN
            YOUT[L]  := YOUT[L]  + TEMP2*PHI[L,I];
            YPOUT[L] := YPOUT[L] + TEMP3*PHI[L,I];
          END;
      END;

    FOR L:=1 TO NEQN DO  YOUT[L] := Y[L] + HI*YOUT[L]; 

  END; (* INTRP *)

