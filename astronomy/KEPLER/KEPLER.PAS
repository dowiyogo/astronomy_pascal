(*-----------------------------------------------------------------------*)
(* KEPLER: calculation of position and velocity for unperturbed          *)
(*         elliptic, parabolic and hyperbolic orbits                     *)
(*                                                                       *)
(*        T0   time of perihelion passage         X,Y,Z     position     *)
(*        T    time                               VX,VY,VZ  velocity     *)
(*        Q    perihelion distance                                       *)
(*        ECC  eccentricity                                              *)
(*        PQR  matrix of Gaussian vectors                                *)
(*        (T0,T in Julian centuries since J2000)                         *)
(*-----------------------------------------------------------------------*)
PROCEDURE KEPLER(T0,T,Q,ECC:REAL;PQR:REAL33;VAR X,Y,Z,VX,VY,VZ:REAL);
  CONST M0=5.0; EPS=0.1;
        KGAUSS = 0.01720209895; DEG = 57.29577951;
  VAR M,DELTA,TAU,INVAX,XX,YY,VVX,VVY: REAL;
  BEGIN
    DELTA := ABS(1.0-ECC);
    INVAX := DELTA / Q;
    TAU   := KGAUSS*36525.0*(T-T0);
    M     := DEG*TAU*SQRT(INVAX*INVAX*INVAX);
    IF ( (M<M0) AND (DELTA<EPS) )
      THEN  PARAB(T0,T,Q,ECC,XX,YY,VVX,VVY)
      ELSE IF (ECC<1.0)
             THEN ELLIP (M,1.0/INVAX,ECC,XX,YY,VVX,VVY)
             ELSE HYPERB(T0,T,1.0/INVAX,ECC,XX,YY,VVX,VVY);
    ORBECL(XX,YY,PQR,X,Y,Z); ORBECL(VVX,VVY,PQR,VX,VY,VZ);
  END;

