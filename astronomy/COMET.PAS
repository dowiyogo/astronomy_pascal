(*-----------------------------------------------------------------------*)
(*                                COMET                                  *)
(*    calculation of unperturbed ephemeris with arbitrary eccentricity   *)
(*                   for comets and minor planets                        *)
(*                          version 93/07/01                             *)
(*-----------------------------------------------------------------------*)

PROGRAM COMET(INPUT,OUTPUT,COMINP);

  uses matlib, pnulib, sphlib, sunlib, keplib, timlib;


  VAR DAY,MONTH,YEAR,NLINE               : INTEGER;
      D,HOUR,T0,Q,ECC,TEQX0,FAC          : REAL;
      MODJD,T,T1,DT,T2,TEQX              : REAL;
      X,Y,Z,VX,VY,VZ,XS,YS,ZS            : REAL;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0 : REAL;
      PQR,A,AS                           : REAL33;
      COMINP                             : TEXT;


(*-----------------------------------------------------------------------*)
(* GETELM: reads orbital elements from file COMINP                       *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETELM (VAR T0,Q,ECC:REAL;VAR PQR:REAL33;VAR TEQX0:REAL);

  VAR INC,LAN,AOP: REAL;

  BEGIN

    WRITELN;
    WRITELN (' COMET: ephemeris calculation for comets and minor planets');
    WRITELN ('                     version 93/07/01                     ');
    WRITELN ('        (C) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    WRITELN;
    WRITELN (' Orbital elements from file COMINP: ');
    WRITELN;

    (* open file for reading *)
    (* RESET(COMINP); *)                              (* Standard Pascal *)
    ASSIGN(COMINP,'COMINP.DAT'); RESET(COMINP);       (* TURBO Pascal    *)

    (* display orbital elements *)
    READLN (COMINP,YEAR,MONTH,D);
    WRITELN('  perihelion time (y m d) ',YEAR:7,MONTH:3,D:6:2);
    READLN(COMINP,Q);
    WRITELN('  perihelion distance (q) ',  Q:14:7,' AU  ');
    READLN(COMINP,ECC);
    WRITELN('  eccentricity (e)        ',ECC:14:7);
    READLN(COMINP,INC);
    WRITELN('  inclination (i)         ',INC:12:5,' deg');
    READLN(COMINP,LAN);
    WRITELN('  long. of ascending node ',LAN:12:5,' deg');
    READLN(COMINP,AOP);
    WRITELN('  argument of perihelion  ',AOP:12:5,' deg');
    READLN(COMINP,TEQX0);
    WRITELN('  equinox                 ',TEQX0:9:2);

    DAY:=TRUNC(D); HOUR:=24.0*(D-DAY);
    T0 := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    TEQX0 := (TEQX0-2000.0)/100.0;
    GAUSVEC(LAN,INC,AOP,PQR);

  END;

(*-----------------------------------------------------------------------*)
(* GETEPH: reads desired period of time and equinox of the ephemeris     *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETEPH(VAR T1,DT,T2,TEQX:REAL);

  VAR YEAR,MONTH,DAY :  INTEGER;
      EQX,HOUR,JD :     REAL;

  BEGIN

    WRITELN;
    WRITELN(' Begin and end of the ephemeris: ');
    WRITELN;
    WRITE  ('  first date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  final date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  step size (dd hh.hh)                      ... ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
    WRITELN;
    WRITE  (' Desired equinox of the ephemeris (yyyy.y)  ... ');
    READLN (EQX);
    TEQX := (EQX-2000.0)/100.0;

    WRITELN; WRITELN;
    WRITELN ('    Date      ET   Sun     l      b     r',
             '        Ra          Dec      Distance ');
    WRITELN (' ':45,'   h  m  s      o  ''  "     (AU) ');

  END;

(*-----------------------------------------------------------------------*)
(* WRTLBR: write L and B in deg,min,sec and R                            *)
(*-----------------------------------------------------------------------*)
PROCEDURE WRTLBR(L,B,R:REAL);
  VAR H,M: INTEGER;
      S  : REAL;
  BEGIN
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); WRITELN(H:5,M:3,TRUNC(S+0.5):3,R:11:6);
  END;
(*------------------------------------------------------------------------*)


BEGIN (* COMET *)


  GETELM (T0,Q,ECC,PQR,TEQX0);  (* read orbital elements              *)
  GETEPH (T1,DT,T2,TEQX);       (* read period of time and equinox    *)

  NLINE := 0;


  PMATECL (TEQX0,TEQX,A);       (* calculate precession matrix         *)

  T := T1;

  REPEAT

    (* date *)

    MODJD := T*36525.0+51544.5;  CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);

    (* ecliptic coordinates of the sun, equinox TEQX        *)

    SUN200 (T,LS,BS,RS);  CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS);  PRECART (AS,XS,YS,ZS);

    (* heliocentric ecliptic coordinates of the comet       *)

    KEPLER (T0,T,Q,ECC,PQR,X,Y,Z,VX,VY,VZ);
    PRECART (A,X,Y,Z);  PRECART (A,VX,VY,VZ);  POLAR (X,Y,Z,R,B,L);

    (* geometric geocentric coordinates of the comet        *)

    X:=X+XS; Y:=Y+YS; Z:=Z+ZS;  DELTA0 := SQRT ( X*X + Y*Y + Z*Z );

    (* first order correction for light travel time         *)

    FAC:=0.00578*DELTA0;  X:=X-FAC*VX;  Y:=Y-FAC*VY;  Z:=Z-FAC*VZ;
    ECLEQU (TEQX,X,Y,Z);  POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

    (* output *)

    WRITE(YEAR:4,'/',MONTH:2,'/',DAY:2,HOUR:6:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);  WRTLBR(RA,DEC,DELTA0);
    NLINE := NLINE+1; IF (NLINE MOD 5) = 0 THEN WRITELN;

    (* next time step                                       *)

    T := T + DT;

  UNTIL (T2<T);


END.
