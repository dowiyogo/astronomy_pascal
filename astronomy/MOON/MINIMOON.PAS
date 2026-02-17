(*-----------------------------------------------------------------------*)
(* MINI_MOON: low precision lunar coordinates (approx. 5'/1')            *)
(*            T  : time in Julian centuries since J2000                  *)
(*                 ( T=(JD-2451545)/36525 )                              *)
(*            RA : right ascension (in h; equinox of date)               *)
(*            DEC: declination (in deg; equinox of date)                 *)
(*-----------------------------------------------------------------------*)
PROCEDURE MINI_MOON (T: REAL; VAR RA,DEC: REAL);
  CONST P2 =6.283185307; ARC=206264.8062;
        COSEPS=0.91748; SINEPS=0.39778;  (* cos/sin(obliquity ecliptic)  *)
  VAR   L0,L,LS,F,D,H,S,N,DL,CB    : REAL;
        L_MOON,B_MOON,V,W,X,Y,Z,RHO: REAL;
  FUNCTION FRAC(X:REAL):REAL; 
    (* with some compilers it may be necessary to replace *)
    (* TRUNC by LONG_TRUNC oder INT if T<-24!             *)
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;
  BEGIN 
    (* mean elements of lunar orbit *)
    L0:=   FRAC(0.606433+1336.855225*T); (* mean longitude Moon (in rev) *)   
    L :=P2*FRAC(0.374897+1325.552410*T); (* mean anomaly of the Moon     *)
    LS:=P2*FRAC(0.993133+  99.997361*T); (* mean anomaly of the Sun      *)
    D :=P2*FRAC(0.827361+1236.853086*T); (* diff. longitude Moon-Sun     *)
    F :=P2*FRAC(0.259086+1342.227825*T); (* mean argument of latitude    *)
    DL := +22640*SIN(L) - 4586*SIN(L-2*D) + 2370*SIN(2*D) +  769*SIN(2*L) 
          -668*SIN(LS)- 412*SIN(2*F) - 212*SIN(2*L-2*D) - 206*SIN(L+LS-2*D)
          +192*SIN(L+2*D) - 165*SIN(LS-2*D) - 125*SIN(D) - 110*SIN(L+LS)
          +148*SIN(L-LS) - 55*SIN(2*F-2*D);
    S := F + (DL+412*SIN(2*F)+541*SIN(LS)) / ARC; 
    H := F-2*D;
    N := -526*SIN(H) + 44*SIN(L+H) - 31*SIN(-L+H) - 23*SIN(LS+H) 
         + 11*SIN(-LS+H) -25*SIN(-2*L+F) + 21*SIN(-L+F);
    L_MOON := P2 * FRAC ( L0 + DL/1296E3 ); (* in rad *)
    B_MOON := ( 18520.0*SIN(S) + N ) / ARC; (* in rad *)
    (* equatorial coordinates *)
    CB:=COS(B_MOON);
    X:=CB*COS(L_MOON); V:=CB*SIN(L_MOON); W:=SIN(B_MOON);
    Y:=COSEPS*V-SINEPS*W; Z:=SINEPS*V+COSEPS*W; RHO:=SQRT(1.0-Z*Z);
    DEC := (360.0/P2)*ARCTAN(Z/RHO); 
    RA  := ( 48.0/P2)*ARCTAN(Y/(X+RHO)); IF RA<0 THEN RA:=RA+24.0;
  END;

