(*-----------------------------------------------------------------------*)
(*                                CONTROL                                *)
(*              confirmation of orbit determination results              *)
(*                            version 93/07/01                           *)
(*-----------------------------------------------------------------------*)

PROGRAM CONTROL(INPUT,OUTPUT,ORBINP,ORBOUT);

  uses matlib, pnulib, sphlib, sunlib, keplib, timlib;

  VAR DAY,MONTH,YEAR,I                       : INTEGER;
      D,MODJD,HOUR,T,T0,Q,ECC,TEQX0,TEQX,FAC : REAL;
      X,Y,Z,VX,VY,VZ,XS,YS,ZS                : REAL;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0     : REAL;
      PQR,A,AS                               : REAL33;
      DATE                                   : ARRAY[1..3] of REAL;
      ORBINP,ORBOUT                          : TEXT;


(*-----------------------------------------------------------------------*)
(* GETELMF: input orbital elements from file ORBOUT                      *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETELMF (VAR T0,Q,ECC:REAL;VAR PQR:REAL33;VAR TEQX0:REAL);

  VAR INC,LAN,AOP: REAL;
      I          : INTEGER;

  BEGIN

    (* RESET(ORBOUT); *)                              (* Standard Pascal *)
    ASSIGN(ORBOUT,'ORBOUT.DAT'); RESET(ORBOUT);       (* Turbo Pascal    *)

    READLN(ORBOUT,YEAR,MONTH,D); DAY:=TRUNC(D); HOUR:=24.0*(D-DAY);
    WRITELN(' Perihelion date (y m d) ',YEAR:4,'/',MONTH:2,'/',D:6:2);
    T0 := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    READLN(ORBOUT,Q);    READLN(ORBOUT,ECC);  READLN(ORBOUT,INC);
    READLN(ORBOUT,LAN);  READLN(ORBOUT,AOP);  READLN(ORBOUT,TEQX0);
    WRITELN(' Perihelion distance (q) ',Q:12:7,' AU');
    WRITELN(' Eccentricity (e)        ',ECC:12:7);
    WRITELN(' Inclination (i)         ',INC:10:5,' deg');
    WRITELN(' Long. of ascending node ',LAN:10:5,' deg');
    WRITELN(' Argument of perihelion  ',AOP:10:5,' deg');
    WRITELN(' Equinox                 ',TEQX0:7:2);

    TEQX0 := (TEQX0-2000)/100;
    GAUSVEC(LAN,INC,AOP,PQR);

  END;

(*-----------------------------------------------------------------------*)
PROCEDURE WRTLBR(L,B,R:REAL);
  VAR H,M: INTEGER;
      S  : REAL;
  BEGIN
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); WRITELN(H:5,M:3,TRUNC(S+0.5):3,R:11:6);
  END;
(*-----------------------------------------------------------------------*)

BEGIN

  (* input                                                  *)

  GETELMF (T0,Q,ECC,PQR,TEQX0); (* orbital elements         *)

  (* RESET(ORBINP); *)                               (* Standard Pascal *)
  ASSIGN(ORBINP,'ORBINP.DAT'); RESET(ORBINP);        (* Turbo Pascal    *)


  READLN(ORBINP);
  FOR I := 1 TO 3 DO
    BEGIN
      READLN(ORBINP,YEAR,MONTH,DAY,HOUR);
      DATE[I] := MJD(DAY,MONTH,YEAR,HOUR);
    END;
  READ(ORBINP,TEQX); TEQX:=(TEQX-2000.0)/100.0;

  WRITELN;
  WRITELN ('    Date      ET   Sun     l      b     r',
           '        Ra          Dec      Distance ');
  WRITELN (' ':45,'   h  m  s      o  ''  "    (AU) ');

  PMATECL (TEQX0,TEQX,A);

  FOR I:=1 TO 3 DO

   BEGIN

    (* date                                                  *)

    MODJD := DATE[I];
    CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);
    T := (MODJD-51544.5)/36525.0;

    (* ecliptic coordinates of the sun, equinox TEQX         *)

    SUN200 (T,LS,BS,RS);
    CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS); PRECART (AS,XS,YS,ZS);

    (* heliocentric ecliptic coordinates of the comet/planet *)

    KEPLER (T0,T,Q,ECC,PQR,X,Y,Z,VX,VY,VZ);
    PRECART (A,X,Y,Z);   PRECART (A,VX,VY,VZ);
    POLAR (X,Y,Z,R,B,L);

    (* geometric geocentric coordinates of the comet         *)

    X:=X+XS; Y:=Y+YS; Z:=Z+ZS;
    DELTA0 := SQRT ( X*X + Y*Y + Z*Z );

    (* first order correction for light travel time          *)

    FAC := 0.00578*DELTA0;
    X := X-FAC*VX;  Y := Y-FAC*VY;  Z := Z-FAC*VZ;
    ECLEQU (TEQX,X,Y,Z);
    POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

    (* output                                               *)

    WRITE(YEAR:5,'/',MONTH:2,'/',DAY:2,HOUR:5:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);
    WRTLBR(RA,DEC,DELTA0);

   END;


END.
