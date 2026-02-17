(*-----------------------------------------------------------------------*)
(* MINI_SUN: low precision solar coordinates (approx. 1')                *)
(*           T  : time in Julian centuries since J2000                   *)
(*                ( T=(JD-2451545)/36525 )                               *)
(*           RA : right ascension (in h; equinox of date)                *)
(*           DEC: declination (in deg; equinox of date)                  *)
(*-----------------------------------------------------------------------*)
PROCEDURE MINI_SUN(T:REAL; VAR RA,DEC: REAL);
  CONST P2  = 6.283185307; COSEPS=0.91748; SINEPS=0.39778;
  VAR   L,M,DL,SL,X,Y,Z,RHO: REAL;
  FUNCTION FRAC(X:REAL):REAL; 
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;
  BEGIN 
    M  := P2*FRAC(0.993133+99.997361*T); 
    DL := 6893.0*SIN(M)+72.0*SIN(2*M);
    L  := P2*FRAC(0.7859453 + M/P2 + (6191.2*T+DL)/1296E3);
    SL := SIN(L);
    X:=COS(L); Y:=COSEPS*SL; Z:=SINEPS*SL; RHO:=SQRT(1.0-Z*Z);
    DEC := (360.0/P2)*ARCTAN(Z/RHO); 
    RA  := ( 48.0/P2)*ARCTAN(Y/(X+RHO)); IF (RA<0) THEN RA:=RA+24.0;
  END;

