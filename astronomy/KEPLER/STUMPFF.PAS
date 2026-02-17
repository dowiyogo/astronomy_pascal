(*-----------------------------------------------------------------------*)
(*  STUMPFF: calculation of Stumpff's functions C1 = sin(E)/E ,          *)
(*           C2 = (1-cos(E))/(E**2) and C3 = (E-sin(E))/(E**3)           *)
(*           for argument E2=E**2                                        *)
(*           (E: eccentric anomaly in radian)                            *)
(*-----------------------------------------------------------------------*)
PROCEDURE STUMPFF(E2:REAL;VAR C1,C2,C3:REAL);
  CONST EPS=1E-12;
  VAR N,ADD: REAL;
  BEGIN
    C1:=0.0; C2:=0.0; C3:=0.0; ADD:=1.0; N:=1.0;
    REPEAT
      C1:=C1+ADD; ADD:=ADD/(2.0*N);
      C2:=C2+ADD; ADD:=ADD/(2.0*N+1.0);
      C3:=C3+ADD; ADD:=-E2*ADD; N:=N+1.0;
    UNTIL ABS(ADD)<EPS;
  END;

