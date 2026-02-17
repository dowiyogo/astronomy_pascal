(*-----------------------------------------------------------------------*)
(* HYPANOM: calculation of the eccentric anomaly H=HYPANOM(MH,ECC) from  *)
(*          mean anomaly MH and eccentricity ECC for                     *)
(*          hyperbolic orbits                                            *)
(*-----------------------------------------------------------------------*)
FUNCTION HYPANOM(MH,ECC:REAL):REAL;
  CONST EPS=1E-10; MAXIT=15;
  VAR   H,F,EXH,SINHH,COSHH: REAL;
        I                  : INTEGER;
  BEGIN
    H:=LN(2.0*ABS(MH)/ECC+1.8); IF (MH<0.0) THEN H:=-H;
    EXH:=EXP(H); SINHH:=0.5*(EXH-1.0/EXH); COSHH:=0.5*(EXH+1.0/EXH);
    F := ECC*SINHH-H-MH; I:=0;
    WHILE ( (ABS(F)>EPS*(1.0+ABS(H+MH))) AND (I<MAXIT) )  DO
      BEGIN
        H := H - F / (ECC*COSHH-1.0);
        EXH:=EXP(H); SINHH:=0.5*(EXH-1.0/EXH); COSHH:=0.5*(EXH+1.0/EXH);
        F := ECC*SINHH-H-MH; I:=I+1;
      END;
    HYPANOM:=H;
    IF (I=MAXIT) THEN  WRITELN(' convergence problems in HYPANOM');
  END;


