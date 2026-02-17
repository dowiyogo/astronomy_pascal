(*-----------------------------------------------------------------------*)
(*                                NEWMOON                                *)
(*           Date of New Moon and ecliptic latitude of the Moon          *)
(*                           version 93/07/01                            *)
(*-----------------------------------------------------------------------*)

PROGRAM NEWMOON (INPUT,OUTPUT);

uses timlib;

CONST D1 = +1236.853086; (* rate of change dD/dT of the mean elongation  *)
                         (* of the Moon from the Sun (in revol./century) *)
      D0 =    +0.827361; (* mean elongation D of the Moon from the sun   *)
                         (* for the epoch J2000 in units of 1rev=360deg  *)

VAR   DAY,MONTH,YEAR,YEAR_CALC : INTEGER;
      HOUR                     : REAL;
      LUNATION_0,LUNATION_I    : INTEGER;
      T_NEW_MOON,MJD_NEW_MOON  : REAL;
      B_MOON                   : REAL;


(*-----------------------------------------------------------------------*)
(* IMPROVE: improves an approximation T for the time of New Moon and     *)
(* finds the ecliptic longitude B of the Moon for that date              *)
(* ( T in julian cent. since J2000, T=(JD-2451545)/36525 )               *)
(*-----------------------------------------------------------------------*)
PROCEDURE IMPROVE ( VAR T,B: REAL);
  CONST P2 =6.283185307;  (* 2*pi *)
        ARC=206264.8062;  (* arcsec per radian *)
  VAR   L,LS,D,F,DLM,DLS,DLAMBDA: REAL;
  (* with some compilers it may be necessary to replace TRUNC *)
  (* by LONG_TRUNC or INT if T<-24!                           *)
  FUNCTION FRAC(X:REAL):REAL;
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;
  BEGIN
    (* mean elements L,LS,D,F of the lunar orbit                         *)
    L  := P2*FRAC(0.374897+1325.552410*T);   (* mean anomaly of the Moon *)
    LS := P2*FRAC(0.993133+  99.997361*T);   (* mean anomaly of the Sun  *)
    D  := P2*(FRAC(0.5+D0+D1*T)-0.5);        (* mean elongation Moon-Sun *)
    F  := P2*FRAC(0.259086+1342.227825*T);   (* long.Moon-long.asc.node  *)
    (* periodic perturbations of the lunar and solar longitude (in")     *)
    DLM := + 22640*SIN(L) - 4586*SIN(L-2*D) + 2370*SIN(2*D) + 769*SIN(2*L)
           - 668*SIN(LS) - 412*SIN(2*F) - 212*SIN(2*L-2*D)
           - 206*SIN(L+LS-2*D) + 192*SIN(L+2*D) - 165*SIN(LS-2*D)
           - 125*SIN(D) - 110*SIN(L+LS) + 148*SIN(L-LS) - 55*SIN(2*F-2*D);
    DLS := + 6893*SIN(LS) + 72*SIN(2*LS);
    (* difference of the true longitudes of Moon and Sun in revolutions  *)
    DLAMBDA  := D / P2   +   ( DLM - DLS) / 1296000.0;
    (* correction for the time of newmmon *)
    T   := T  - DLAMBDA / D1;
    (* ecliptic latitude B of the moon (in deg) *)
    B  := ( + 18520.0*SIN(F+DLM/ARC) - 526*SIN(F-2*D) ) / 3600.0;
  END;
(*-----------------------------------------------------------------------*)

BEGIN (* main program *)

  WRITELN;
  WRITELN (' NEWMOON: Date of New Moon and ecliptic latitude of the Moon');
  WRITELN ('                       version 93/07/01                     ');
  WRITELN ('         (c) 1993 Thomas Pfleger, Oliver Montenbruck        ');
  WRITELN;
  WRITE   (' Dates of New Moon for the year '); READLN(YEAR_CALC);
  WRITELN;
  WRITELN (' Date ':16,'UT':7,'Latitude':12); WRITELN;
  WRITELN (' h':23,'o':9);

  LUNATION_0 := TRUNC( D1 * (YEAR_CALC-2000)/100 );
  FOR LUNATION_I := LUNATION_0 TO LUNATION_0 + 13 DO
    BEGIN
      T_NEW_MOON   := ( LUNATION_I - D0 ) / D1;
      IMPROVE ( T_NEW_MOON, B_MOON );
      IMPROVE ( T_NEW_MOON, B_MOON );
      MJD_NEW_MOON := 36525.0*T_NEW_MOON + 51544.5;
      CALDAT ( MJD_NEW_MOON, DAY,MONTH,YEAR,HOUR );
      IF YEAR=YEAR_CALC THEN
        WRITELN(YEAR:12,'/',MONTH:2,'/',DAY:2,HOUR:6:1,B_MOON:9:1)
    END;
  WRITELN;

END.
