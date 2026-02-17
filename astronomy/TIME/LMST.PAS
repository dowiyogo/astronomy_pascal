(*-----------------------------------------------------------------------*)
(* LMST: local mean sidereal time                                        *)
(*-----------------------------------------------------------------------*)
FUNCTION LMST(MJD,LAMBDA:REAL):REAL;
  VAR MJD0,T,UT,GMST: REAL;
  FUNCTION FRAC(X:REAL):REAL;
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;
  BEGIN
    (* MJD0:=TRUNC(MJD); *)       (* Standard Pascal  *)
    MJD0:=INT(MJD);               (* TURBO Pascal     *)
    (* MJD0:=LONG_TRUNC(MJD); *)  (* ST Pascal plus   *)
    UT:=(MJD-MJD0)*24; T:=(MJD0-51544.5)/36525.0;
    GMST:=6.697374558 + 1.0027379093*UT
            +(8640184.812866+(0.093104-6.2E-6*T)*T)*T/3600.0;
    LMST:=24.0*FRAC( (GMST-LAMBDA/15.0) / 24.0 );
  END;

