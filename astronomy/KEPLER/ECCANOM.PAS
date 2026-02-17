(*-----------------------------------------------------------------------*)
(* ECCANOM: calculation of the eccentric anomaly E=ECCANOM(MAN,ECC)      *)
(*          from the mean anomaly MAN and the eccentricity ECC.          *)
(*          (solution of Kepler's equation by Newton's method)           *)
(*          (E, MAN in degrees)                                          *)
(*-----------------------------------------------------------------------*)
FUNCTION ECCANOM(MAN,ECC:REAL):REAL;
  CONST PI=3.141592654; TWOPI=6.283185308; RAD=0.0174532925199433;
        EPS = 1E-11; MAXIT = 15;
  VAR M,E,F: REAL;
      I    : INTEGER;
  BEGIN
    M:=MAN/360.0;  M:=TWOPI*(M-TRUNC(M)); IF M<0 THEN M:=M+TWOPI;
    IF (ECC<0.8) THEN E:=M ELSE E:=PI;
    F := E - ECC*SIN(E) - M; I:=0;
    WHILE ( (ABS(F)>EPS) AND (I<MAXIT) ) DO
      BEGIN
        E := E - F / (1.0-ECC*COS(E));  F := E-ECC*SIN(E)-M; I:=I+1;
      END;
    ECCANOM:=E/RAD;
    IF (I=MAXIT) THEN  WRITELN(' convergence problems in ECCANOM');
  END;

