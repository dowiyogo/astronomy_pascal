(*-----------------------------------------------------------------------*)
(* PARAB: calculation of position and velocity for                       *)
(*        parabolic and near parabolic orbits according to Stumpff       *)
(*                                                                       *)
(*        T0   time of perihelion passage         X,Y    position        *)
(*        T    time                              VX,VY  velocity         *)
(*        Q    perihelion distance                                       *)
(*        ECC  eccentricity                                              *)
(*        (T0,T in julian centuries since J2000)                         *)
(*-----------------------------------------------------------------------*)
PROCEDURE PARAB(T0,T,Q,ECC:REAL;VAR X,Y,VX,VY:REAL);
  CONST EPS    = 1E-9;
        KGAUSS = 0.01720209895;
        MAXIT  = 15;
  VAR E2,E20,FAC,C1,C2,C3,K,TAU,A,U,U2: REAL;
      R: REAL;
      I: INTEGER;
  BEGIN
    E2:=0.0;  FAC:=0.5*ECC;  I:=0;
    K   := KGAUSS / SQRT(Q*(1.0+ECC));
    TAU := KGAUSS * 36525.0*(T-T0);
    REPEAT
      I:=I+1;
      E20:=E2;
      A:=1.5*SQRT(FAC/(Q*Q*Q))*TAU;  A:=CUBR(SQRT(A*A+1.0)+A);
      U:=A-1.0/A;  U2:=U*U;  E2:=U2*(1.0-ECC)/FAC;
      STUMPFF(E2,C1,C2,C3); FAC:=3.0*ECC*C3;
    UNTIL (ABS(E2-E20)<EPS)OR(I>MAXIT);
    IF (I=MAXIT) THEN  WRITELN(' convergence problems in PARAB');
    R :=Q*(1.0+U2*C2*ECC/FAC);
    X :=Q*(1.0-U2*C2/FAC);          VY:= K*(X/R+ECC);
    Y :=Q*SQRT((1.0+ECC)/FAC)*U*C1; VX:=-K*Y/R;
  END;

