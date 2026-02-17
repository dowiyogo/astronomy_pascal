(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* STEP                                                                  *)
(*                                                                       *)
(* Basic integration routine for the variable order and stepsize code    *)
(* of Shampine and Gordon.                                               *)
(*                                                                       *)
(*                                                                       *)
(*   X                                                                   *)
(*   Y                                                                   *)
(*   F       Function to be integrated                                   *)
(*   NEQN    Number of differential equations                            *)
(*   H                                                                   *)
(*   EPS                                                                 *)
(*   WT                                                                  *)
(*   START                                                               *)
(*   HOLD                                                                *)
(*   K                                                                   *)
(*   KOLD                                                                *)
(*   CRASH   TRUE if STEP was aborted (tolerances too small)             *)
(*   PHI                                                                 *)
(*   P                                                                   *)
(*   YP                                                                  *)
(*   PSI                                                                 *)
(*   ALPHA                                                               *)
(*   BETA                                                                *)
(*   SIG                                                                 *)
(*   V                                                                   *)
(*   W                                                                   *)
(*   G                                                                   *)
(*   PHASE1                                                              *)
(*   NS                                                                  *)
(*   NORND                                                               *)
(*                                                                       *)
(*                                                                       *)
(* Note:                                                                 *)
(*                                                                       *)
(* Replace "F:DE_FUNC;" by                                               *)
(*   "PROCEDURE F(X:REAL;Y:DE_EQN_VECTOR;VAR YP:DE_EQN_VECTOR);"         *)
(* for use with standard Pascal.                                         *)
(*                                                                       *)
(* A detailed description of the method is given in                      *)
(*                                                                       *)
(*    L. F. Shampine, M. K. Gordon; Computer Solution of ordinary        *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)


PROCEDURE STEP ( VAR X              :  REAL; 
                 VAR Y              :  DE_EQN_VECTOR;
                 F                  :  DE_FUNC;
                 VAR NEQN           :  INTEGER;
                 VAR H,EPS          :  REAL;
                 VAR WT             :  DE_EQN_VECTOR;
                 VAR START          :  BOOLEAN;
                 VAR HOLD           :  REAL;
                 VAR K,KOLD         :  INTEGER;
                 VAR CRASH          :  BOOLEAN;
                 VAR PHI            :  DE_PHI_ARRAY;
                 VAR P,YP           :  DE_EQN_VECTOR;
                 VAR PSI,ALPHA,BETA :  DE_12_VECTOR;
                 VAR SIG            :  DE_13_VECTOR;
                 VAR V,W            :  DE_12_VECTOR;
                 VAR G              :  DE_13_VECTOR;
                 VAR PHASE1         :  BOOLEAN;
                 VAR NS             :  INTEGER;
                 VAR NORND          :  BOOLEAN            );


  LABEL  999;   (* Error exit at end of procedure *)


  VAR I, IFAIL, IM1, IP1, IQ, J, KM1, KM2, KNEW, KP1, KP2  :  INTEGER;
      L, LIMIT1, LIMIT2, NSM2, NSP1, NSP2                  :  INTEGER;
      ABSH, ERK, ERKM1, ERKM2, ERKP1, ERR, HNEW            :  REAL;
      P5EPS, R, REALI, REALNS, RHO, ROUND, SUM, TAU, TEMP1 :  REAL;
      TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, XOLD              :  REAL;
      TWO, GSTR                                            :  DE_13_VECTOR;
      SUCCESS                                              :  BOOLEAN;

  (* SIGN: computes sign(B)*abs(A) *)

  FUNCTION SIGN(A,B: REAL): REAL; 
    BEGIN IF (B>=0.0) THEN SIGN:=ABS(A) ELSE SIGN:=-ABS(A); END;

  (* MIN: computes the minimum of A and B *)

  FUNCTION MIN(A,B: REAL): REAL;
    BEGIN IF A<B THEN MIN:=A ELSE MIN:=B; END;

  (* MAX: computes the maximum of A and B *)

  FUNCTION MAX(A,B: REAL): REAL;
    BEGIN IF A>B THEN MAX:=A ELSE MAX:=B; END;

  BEGIN (* STEP *)


    (* explicit initialization of arrays TWO and GSTR *)

    TWO[1]:=2.0;
    FOR I:=2 TO 13 DO TWO[I]:=2.0*TWO[I-1];

    GSTR[ 1]:=0.5;     GSTR[ 2]:=0.0833;  GSTR[ 3]:=0.0417;
    GSTR[ 4]:=0.0264;  GSTR[ 5]:=0.0188;  GSTR[ 6]:=0.0143;
    GSTR[ 7]:=0.0114;  GSTR[ 8]:=0.00936; GSTR[ 9]:=0.00789;
    GSTR[10]:=0.00679; GSTR[11]:=0.00592; GSTR[12]:=0.00524;
    GSTR[13]:=0.00468; 


    (*                                                                   *)
    (* Begin block 0                                                     *)
    (*                                                                   *)
    (* Check if step size or error tolerance is too small for machine    *)
    (* precision.  If first step, initialize PHI array and estimate a    *)
    (* starting step size. If step size is too small, determine an       *)
    (* acceptable one.                                                   *)
    (*                                                                   *)

    IF ( ABS(H) < FOURU*ABS(X) ) THEN
      BEGIN
        H := SIGN(FOURU*ABS(X),H);
        CRASH := TRUE;
        GOTO 999;   (* Exit *)
      END; 

    P5EPS := 0.5*EPS;
    CRASH := FALSE;
    G[1]  := 1.0;
    G[2]  := 0.5;
    SIG[1] := 1.0;

    IFAIL := 0;

    (* If error tolerance is too small, increase it to an *)
    (* acceptable value.                                  *)

    ROUND := 0.0;
    FOR L:=1 TO NEQN DO  ROUND := ROUND + (Y[L]*Y[L])/(WT[L]*WT[L]);
    ROUND := TWOU*SQRT(ROUND);
    IF P5EPS<ROUND THEN
      BEGIN
        EPS := 2.0*ROUND*(1.0+FOURU);
        CRASH := TRUE;
        GOTO 999;
      END; 


    IF (START) THEN 

      BEGIN

        (* Initialize. Compute appropriate step size for first step. *)
 
        F(X,Y,YP);
        SUM := 0.0;
        FOR L:=1 TO NEQN DO 
          BEGIN 
            PHI[L,1] := YP[L];
            PHI[L,2] := 0.0;
            SUM := SUM + (YP[L]*YP[L])/(WT[L]*WT[L]);
          END; 
        SUM := SQRT(SUM);
        ABSH := ABS(H);
        IF (EPS<16.0*SUM*H*H) THEN ABSH:=0.25*SQRT(EPS/SUM);
        H := SIGN(MAX(ABSH,FOURU*ABS(X)),H);
        HOLD := 0.0;
        K := 1;
        KOLD := 0;
        START := FALSE;
        PHASE1 := TRUE;
        NORND := TRUE;
        IF (P5EPS<=100.0*ROUND) THEN
          BEGIN
            NORND := FALSE;
            FOR L:=1 TO NEQN DO PHI[L,15]:=0.0;
          END;

      END;

    (*                                                                   *)
    (* End block 0                                                       *)
    (*                                                                   *)


    (*                                                                   *)
    (* Repeat blocks 1, 2 (and 3) until step is successful               *)
    (*                                                                   *)

    REPEAT

      (*                                                                 *)
      (* Begin block 1                                                   *)
      (*                                                                 *)
      (* Compute coefficients of formulas for this step. Avoid computing *)
      (* those quantities not changed when step size is not changed.     *)
      (*                                                                 *)
  
      KP1 := K+1;
      KP2 := K+2;
      KM1 := K-1;
      KM2 := K-2;
  
      (* NS is the number of steps taken with size H, including the *)
      (* current one.  When K<NS, no coefficients change.           *)
  
      IF (H<>HOLD) THEN NS:=0;
      IF (NS<=KOLD) THEN NS:=NS+1;
      NSP1 := NS+1;
  
  
      IF (K>=NS) THEN 
  
        BEGIN
  
          (* Compute those components of ALPHA[*],BETA[*],PSI[*],SIG[*] *)
          (* which are changed                                          *)
  
          BETA[NS] := 1.0;
          REALNS := NS;
          ALPHA[NS] := 1.0/REALNS;
          TEMP1 := H*REALNS;
          SIG[NSP1] := 1.0;
          IF (K>=NSP1) THEN
            FOR I:=NSP1 TO K DO 
              BEGIN
                IM1 := I-1;
                TEMP2 := PSI[IM1];
                PSI[IM1] := TEMP1;
                BETA[I] := BETA[IM1]*PSI[IM1]/TEMP2;
                TEMP1 := TEMP2 + H;
                ALPHA[I] := H/TEMP1;
                REALI := I;
                SIG[I+1] := REALI*ALPHA[I]*SIG[I];
              END;
          PSI[K] := TEMP1;
  
          (* Compute coefficients G[*]; initialize V[*] and set W[*]. *)
  
          IF (NS>1) 
  
            THEN
  
              BEGIN

                (* If order was raised, update diagonal part of V[*] *)

                IF (K>KOLD) THEN 
                  BEGIN
                    TEMP4 := K*KP1;
                    V[K] := 1.0/TEMP4;
                    NSM2 := NS-2;
                    FOR J:=1 TO NSM2 DO 
                      BEGIN
                        I := K-J;
                        V[I] := V[I] - ALPHA[J+1]*V[I+1];
                      END;
                  END;

                (* Update V[*] and set W[*] *)

                LIMIT1 := KP1 - NS;
                TEMP5 := ALPHA[NS];
                FOR IQ:=1 TO LIMIT1 DO
                  BEGIN
                    V[IQ] := V[IQ] - TEMP5*V[IQ+1];
                    W[IQ] := V[IQ];
                  END;
                G[NSP1] := W[1];

              END
  
            ELSE 
    
              FOR IQ:=1 TO K DO 
                BEGIN
                  TEMP3 := IQ*(IQ+1);
                  V[IQ] := 1.0/TEMP3;
                  W[IQ] := V[IQ];
                END;
  
  
          (* Compute the G[*] in the work vector W[*] *)
  
          NSP2 := NS + 2;
          IF (KP1>=NSP2) THEN 
            FOR I:=NSP2 TO KP1 DO 
              BEGIN
                LIMIT2 := KP2 - I;
                TEMP6 := ALPHA[I-1];
                FOR IQ:=1 TO LIMIT2 DO  W[IQ] := W[IQ] - TEMP6*W[IQ+1];
                G[I] := W[1];
              END;
  
  
        END; (* if K>=NS *) 
  
      (*                                                                 *)
      (* End block 1                                                     *)
      (*                                                                 *)
  
  
  
      (*                                                                 *)
      (* Begin block 2                                                   *)
      (*                                                                 *)
      (* Predict a solution P[*], evaluate derivatives using predicted   *)
      (* solution, estimate local error at order K and errors at orders  *)
      (* K, K-1, K-2 as if constant step size were used.                 *)
      (*                                                                 *)
  
      (* Change PHI to PHI star *)
  
      IF (K>=NSP1) THEN
        FOR I:=NSP1 TO K DO 
          BEGIN
            TEMP1 := BETA[I];
            FOR L:=1 TO NEQN DO PHI[L,I] := TEMP1*PHI[L,I];
          END;
  
      (* Predict solution and differences *)
  
      FOR L:=1 TO NEQN DO 
        BEGIN
          PHI[L,KP2] := PHI[L,KP1];
          PHI[L,KP1] := 0.0;
          P[L] := 0.0;
        END;
      FOR J:=1 TO K DO 
        BEGIN 
          I := KP1 - J;
          IP1 := I+1;
          TEMP2 := G[I];
          FOR L:=1 TO NEQN DO 
            BEGIN
              P[L] := P[L] + TEMP2*PHI[L,I];
              PHI[L,I] := PHI[L,I] + PHI[L,IP1];
            END;
        END;
      IF (NORND) 
        THEN 
          FOR L:=1 TO NEQN DO  P[L] := Y[L] + H*P[L]
        ELSE
          FOR L:=1 TO NEQN DO
            BEGIN
              TAU := H*P[L] - PHI[L,15];
              P[L] := Y[L] + TAU;
              PHI[L,16] := (P[L] - Y[L]) - TAU;
            END;
      XOLD := X;
      X := X + H;
      ABSH := ABS(H);
      F(X,P,YP);
  
      (* Estimate errors at orders K,K-1,K-2 *)
  
      ERKM2 := 0.0;
      ERKM1 := 0.0;
      ERK := 0.0;
      FOR L:=1 TO NEQN DO 
        BEGIN
          TEMP3 := 1.0/WT[L];
          TEMP4 := YP[L] - PHI[L,1];
          IF (KM2>0) THEN 
             ERKM2 := ERKM2 + ((PHI[L,KM1]+TEMP4)*TEMP3)
                             *((PHI[L,KM1]+TEMP4)*TEMP3);
          IF (KM2>=0) THEN 
             ERKM1 := ERKM1 + ((PHI[L,K]+TEMP4)*TEMP3)
                             *((PHI[L,K]+TEMP4)*TEMP3);
          ERK := ERK + (TEMP4*TEMP3)*(TEMP4*TEMP3);
        END;
      IF (KM2>0) THEN ERKM2 := ABSH*SIG[KM1]*GSTR[KM2]*SQRT(ERKM2);
      IF (KM2>=0) THEN ERKM1 := ABSH*SIG[K]*GSTR[KM1]*SQRT(ERKM1);
      TEMP5 := ABSH*SQRT(ERK);
      ERR := TEMP5*(G[K]-G[KP1]);
      ERK := TEMP5*SIG[KP1]*GSTR[K];
      KNEW := K;
  
      (* Test if order should be lowered *)
  
      IF (KM2>0) THEN 
        IF (MAX(ERKM1,ERKM2)<=ERK) THEN  KNEW:=KM1;
      IF (KM2=0) THEN 
        IF (ERKM1<=0.5*ERK) THEN KNEW:=KM1;
  
      (*                                                                 *)
      (* End block 2                                                     *)
      (*                                                                 *)
  

  
      (*                                                                 *)
      (* If step is successful continue with block 4, otherwise repeat   *)
      (* blocks 1 and 2 after executing block 3                          *)
      (*                                                                 *)
  
      SUCCESS := (ERR<=EPS);
  
  
      IF (NOT SUCCESS) THEN 
  
        BEGIN 
  
          (*                                                             *)
          (* Begin block 3                                               *)
          (*                                                             *)
          (* The step is unsuccessful. Restore X, PHI[*,*], PSI[*]. If   *)
          (* 3rd consecutive failure, set order to 1. If step fails more *)
          (* than 3 times, consider an optimal step size. Double error   *)
          (* tolerance and return if estimated step size is too small    *)
          (* for machine precision.                                      *)
          (*                                                             *)
      
          (* Restore X, PHI[*,*] and PSI[*] *)
      
          PHASE1 := FALSE; 
          X := XOLD;
          FOR I:=1 TO K DO 
            BEGIN 
              TEMP1 := 1.0/BETA[I];
              IP1 := I+1;
              FOR L:=1 TO NEQN DO  PHI[L,I]:=TEMP1*(PHI[L,I]-PHI[L,IP1]);
            END;
          IF (K>=2) THEN 
            FOR I:=2 TO K DO  PSI[I-1] := PSI[I] - H;
      
          (* On third failure, set order to one. *)
          (* Thereafter, use optimal step size   *)
      
          IFAIL := IFAIL + 1;
          TEMP2 := 0.5;
          IF (IFAIL>3) THEN 
            IF (P5EPS < 0.25*ERK) THEN TEMP2 := SQRT(P5EPS/ERK);
          IF (IFAIL>=3) THEN  KNEW := 1;
          H := TEMP2*H;
          K := KNEW;
          IF (ABS(H)<FOURU*ABS(X)) THEN 
            BEGIN
              CRASH := TRUE;
              H := SIGN(FOURU*ABS(X),H);
              EPS := EPS + EPS;
              GOTO 999;   (* Exit *)
            END;
      
          (*                                                             *)
          (* End block 3, return to start of block 1                     *)
          (*                                                             *)
  
        END;  (* end if(successful) *)
  

    UNTIL (SUCCESS);



    (*                                                                   *)
    (* Begin block 4                                                     *)
    (*                                                                   *)
    (* The step is successful. Correct the predicted solution, evaluate  *)
    (* the derivatives using the corrected solution and update the       *)
    (* differences. Determine best order and step size for next step.    *)
    (*                                                                   *)

    KOLD := K;
    HOLD := H;

    (* Correct and evaluate *)

    TEMP1 := H*G[KP1];
    IF (NORND) 
      THEN 
        FOR L:=1 TO NEQN DO  Y[L] := P[L] + TEMP1*(YP[L] - PHI[L,1])
      ELSE 
        FOR L:=1 TO NEQN DO
          BEGIN
            RHO := TEMP1*(YP[L] - PHI[L,1]) - PHI[L,16];
            Y[L] := P[L] + RHO;
            PHI[L,15] := (Y[L] - P[L]) - RHO;
          END;
    F(X,Y,YP);

    (* Update differences for next step *)

    FOR L:=1 TO NEQN DO
      BEGIN
        PHI[L,KP1] := YP[L] - PHI[L,1];
        PHI[L,KP2] := PHI[L,KP1] - PHI[L,KP2];
      END;
    FOR I:=1 TO K DO
      FOR L:=1 TO NEQN DO
        PHI[L,I] := PHI[L,I] + PHI[L,KP1];

    (* Estimate error at order K+1 unless               *)
    (* - in first phase when always raise order,        *)
    (* - already decided to lower order,                *)
    (* - step size not constant so estimate unreliable  *)

    ERKP1 := 0.0;
    IF ( (KNEW=KM1) OR (K=12) ) THEN  PHASE1:=FALSE;

    IF ( PHASE1 ) 

      THEN 

        BEGIN
          K := KP1;
          ERK := ERKP1;
        END

      ELSE

        BEGIN

          IF ( KNEW=KM1 ) 

            THEN  (* lower order *)

              BEGIN  
                K := KM1;
                ERK := ERKM1;
              END

            ELSE

              IF ( KP1<=NS ) THEN 
                BEGIN
                  FOR L:=1 TO NEQN DO 
                    ERKP1 := ERKP1 + (PHI[L,KP2]/WT[L])*(PHI[L,KP2]/WT[L]);
                  ERKP1 := ABSH*GSTR[KP1]*SQRT(ERKP1);
                  (* Using estimated error at order K+1, determine *)
                  (* appropriate order for next step               *)
                  IF (K>1) 
                    THEN
                      IF ( ERKM1<=MIN(ERK,ERKP1)) 
                        THEN  (* lower order *)
                          BEGIN K:=KM1; ERK:=ERKM1 END
                        ELSE 
                          BEGIN
                            IF ( (ERKP1<ERK) AND (K<>12) ) THEN 
                               (* raise order *)
                               BEGIN K:=KP1; ERK:=ERKP1; END;
                          END
                    ELSE
                      IF ( ERKP1<0.5*ERK) THEN  (* raise order *)
                        BEGIN
                          (* Here ERKP1 < ERK < max(ERKM1,ERKM2) else    *)
                          (* order would have been lowered in block 2.   *)
                          (* Thus order is to be raised                  *)
                          K := KP1;
                          ERK := ERKP1;
                        END;
                  (* end if  K>1 *)
                END; (* if KP1<=NS *)

             (* end if KNEW=KM1 *)

        END;  (* if PHASE1 *)
 

    (* With new order determine appropriate step size for next step *)

    IF ( (PHASE1) OR (P5EPS>=ERK*TWO[K+1]) ) 
      THEN
        HNEW := H + H
      ELSE 
        BEGIN
          IF ( P5EPS<ERK ) 
            THEN 
              BEGIN
                TEMP2 := K+1;
                R := EXP ( LN(P5EPS/ERK) * (1.0/TEMP2) );
                HNEW := ABSH*MAX(0.5,MIN(0.9,R));
                HNEW := SIGN(MAX(HNEW,FOURU*ABS(X)),H);
              END
            ELSE
              HNEW := H;
        END;

    H := HNEW;


    (*                                                                   *)
    (* End block 4                                                       *)
    (*                                                                   *)


  999: 

  END;  (* STEP *)

