(*-----------------------------------------------------------------------*)
(* FIND_ETA: determines the sector/triangle ratio                        *)
(*           from two positions and the time difference                  *)
(*-----------------------------------------------------------------------*)
FUNCTION FIND_ETA ( RA,RB: VECTOR; TAU: REAL ): REAL;

  CONST DELTA = 1.0E-9;  MAXIT = 30;

  VAR   KAPPA,M,L,SA,SB,ETA_MIN,ETA1,ETA2,F1,F2,D_ETA: REAL;
        I: INTEGER;

  (* F(eta) = 1 - eta + (m/eta**2)*W(m/eta**2-l) *)
  FUNCTION F ( ETA,M,L: REAL ): REAL;
    CONST EPS =1.0E-10;
    VAR W,WW,A,S,N,G,E: REAL;
    BEGIN
      W := M/(ETA*ETA)-L;
      IF (ABS(W)<0.1)
        THEN  (* series expansion *)
          BEGIN
            A:=4.0/3.0; WW:=A; N:=0.0;
            REPEAT
              N:=N+1; A:=A*W*(N+2.0)/(N+1.5); WW:=WW+A;
            UNTIL ABS(A)<EPS;
          END
        ELSE
          IF (W>0)
            THEN (* W=(2g-sin2g)/(sin(g)**3), g=2*arcsin(sqrt(w)) *)
              BEGIN
                G := 2.0*ARCTAN(SQRT(W/(1.0-W)));  S := SIN(G);
                WW := (2.0*G-SIN(2.0*G))/(S*S*S);
              END
            ELSE (* W=(sinh2g-2g)/(sinh(g)**3), g=2*arsinh(sqrt(-w)) *)
              BEGIN
                G := 2.0*LN(SQRT(-W)+SQRT(1.0-W));
                E := EXP(G); S:=0.5*(E-1.0/E); E:=E*E;
                WW := (0.5*(E-1.0/E)-2.0*G)/(S*S*S);
              END;
       F := 1.0-ETA+(W+L)*WW;
     END; (* FIND_ETA.F *)

  BEGIN

    SA := NORM(RA);  SB := NORM(RB);  KAPPA := SQRT(2.0*(SA*SB+DOT(RA,RB)));
    M := TAU*TAU / (KAPPA*KAPPA*KAPPA);   L := (SA+SB)/(2.0*KAPPA) - 0.5;
    ETA_MIN := SQRT(M/(L+1.0));

    (* start with Hansen's approximation *)
    ETA2 := ( 12.0 + 10.0*SQRT(1.0+(44.0/9.0)*M/(L+5.0/6.0)) ) / 22.0;
    ETA1 := ETA2 + 0.1;   F1 := F(ETA1,M,L);   F2 := F(ETA2,M,L);  I := 0;

    (* secant method *)
    WHILE ( (ABS(F2-F1)>DELTA) AND (I<MAXIT) ) DO
      BEGIN
        D_ETA:=-F2*(ETA2-ETA1)/(F2-F1);  ETA1:=ETA2; F1:=F2;
        WHILE (ETA2+D_ETA<=ETA_MIN) DO D_ETA:=0.5*D_ETA;
        ETA2:=ETA2+D_ETA;  F2:=F(ETA2,M,L);  I:=I+1;
      END;
    IF (I=MAXIT) THEN  WRITELN(' convergence problems in FIND_ETA');
    FIND_ETA := ETA2;

  END;


