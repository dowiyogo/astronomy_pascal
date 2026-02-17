(*-----------------------------------------------------------------------*)
(*                                  COCO                                 *)
(*                         coordinate conversion                         *)
(*                            version 93/07/01                           *)
(*-----------------------------------------------------------------------*)

PROGRAM COCO(INPUT,OUTPUT);

  uses matlib, pnulib, sphlib, sunlib, timlib;

  VAR X,Y,Z,XS,YS,ZS: REAL;
      T,TEQX,TEQXN  : REAL;
      LS,BS,RS      : REAL;
      A             : REAL33;
      ECLIPT        : BOOLEAN;
      MODE          : CHAR;


(*-----------------------------------------------------------------------*)

PROCEDURE GETEQX(VAR TEQX:REAL);
  BEGIN
    WRITE  (' equinox (yyyy) ? ');
    READLN (TEQX); TEQX := (TEQX-2000.0)/100.0;
  END;

(*-----------------------------------------------------------------------*)

PROCEDURE GETDAT(VAR T:REAL);
  VAR D,M,Y  : INTEGER;
      HOUR,JD: REAL;
  BEGIN
    WRITE (' Date (year month day hour) ?      ');
    READLN (Y,M,D,HOUR);
    JD := MJD(D,M,Y,HOUR) + 2400000.5;
    WRITELN; WRITELN (' JD',JD:13:4); WRITELN;
    T := (JD-2451545.0) / 36525.0;
  END;

(*-----------------------------------------------------------------------*)

PROCEDURE GETINP (VAR X,Y,Z,TEQX:REAL;VAR ECLIPT:BOOLEAN);

  VAR I,D,M  : INTEGER;
      L,B,R,S: REAL;

  BEGIN (* GETINP *)

    WRITELN;
    WRITELN('                COCO: coordinate conversion           ');
    WRITELN('                      version 93/07/01                ');
    WRITELN('        (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    WRITELN;
    WRITELN (' Coordinate input: please select format required ');
    WRITELN;
    WRITELN ('  1  ecliptic   cartesian     2  ecliptic   polar');
    WRITELN ('  3  equatorial cartesian     4  equatorial polar');
    WRITELN;
    WRITE   ('  '); READLN (I); WRITELN;

    CASE I OF
      1: BEGIN
           WRITE (' Coordinates (x y z) ?  '); READLN(X,Y,Z); ECLIPT:=TRUE;
         END;
      2: BEGIN
           WRITE (' Coordinates (L (o '' ")  B (o '' ")  R) ?  ');
           READ(D,M,S); DDD(D,M,S,L); READLN(D,M,S,R); DDD(D,M,S,B);
           CART(R,B,L,X,Y,Z); ECLIPT:=TRUE;
         END;
      3: BEGIN
           WRITE (' Coordinates (x y z) ?  '); READLN(X,Y,Z);ECLIPT:=FALSE;
         END;
      4: BEGIN
           WRITE (' Coordinates (RA (h m s)  DEC (o '' ")  R) ?  ');
           READ(D,M,S); DDD(D,M,S,L); READLN(D,M,S,R); DDD(D,M,S,B);
           L:=L*15.0; CART(R,B,L,X,Y,Z); ECLIPT:=FALSE;
         END;
      END; (* CASE *)

    GETEQX (TEQX); (* read equinox *)

  END; (* GETINP *)

(*-----------------------------------------------------------------------*)

PROCEDURE RESULT (X,Y,Z: REAL; ECLIPT: BOOLEAN);

  VAR L,B,R,S: REAL;
      D,M    : INTEGER;

  BEGIN (* RESULT *)

    WRITELN; WRITELN (' (x,y,z) = (',x:13:8,',',y:13:8,',',z:13:8,')');
    WRITELN;

    POLAR (X,Y,Z,R,B,L);
    IF ECLIPT
      THEN
        BEGIN
          WRITELN (' ':5,'   o  ''  " ',' ':8,'   o  ''  " ');
          DMS(L,D,M,S); WRITE(' L = ',D:3,M:3,S:5:1,' ':3);
          DMS(B,D,M,S); WRITE(' B = ',D:3,M:3,S:5:1,' ':3);
        END
      ELSE
        BEGIN
          WRITELN (' ':5,'   h  m  s ',' ':10,'   o  ''  " ');
          DMS(L/15,D,M,S); WRITE(' RA = ',D:2,M:3,S:5:1,' ':3);
          DMS(B,D,M,S);    WRITE(' DEC = ',D:3,M:3,S:5:1,' ':3);
        END;

    WRITELN (' R = ',R:12:8); WRITELN; WRITELN;

  END; (* RESULT *)

(*-----------------------------------------------------------------------*)

BEGIN (* COCO *)

  GETINP (X,Y,Z,TEQX,ECLIPT);
  RESULT (X,Y,Z,ECLIPT);

  REPEAT

    WRITE (' Command (?=Help): ');  READLN (MODE);   WRITELN;

    IF MODE IN ['?','P','p','E','e','H','h','G','g'] THEN
    CASE MODE OF

      '?'    : BEGIN (* help *)
                 WRITELN;
                 WRITELN ('   E: equatorial <-> ecliptic                ');
                 WRITELN ('   P: -> precession       G: -> geocentric   ');
                 WRITELN ('   H: -> heliocentric     S: -> STOP         ');
                 WRITELN;
               END;

      'P','p': BEGIN (* precession *)
                 WRITE (' New');
                 GETEQX(TEQXN); (* read new equinox *)
                 IF ECLIPT THEN PMATECL(TEQX,TEQXN,A)
                           ELSE PMATEQU(TEQX,TEQXN,A);
                 PRECART(A,X,Y,Z);
                 TEQX := TEQXN;
                 WRITELN;
                 WRITELN (' Coordinates referred to equinox T =',
                          TEQX:13:10);
               END;

      'E','e': BEGIN (* ecliptic <-> equatorial *)
                 WRITELN;
                 IF (ECLIPT) THEN
                   BEGIN ECLEQU(TEQX,X,Y,Z); WRITE(' Equatorial'); END
                 ELSE
                   BEGIN EQUECL(TEQX,X,Y,Z); WRITE(' Ecliptic'); END;
                 WRITELN (' coordinates: ');
                 ECLIPT := NOT ECLIPT;
               END;

      'G','g', (* -> geocentric coordinates   *)
      'H','h': (* -> heliocentric coordinates *)
               BEGIN
                 GETDAT(T); (* read date *)
                 SUN200(T,LS,BS,RS);
                 CART(RS,BS,LS,XS,YS,ZS);
                 PMATECL(T,TEQX,A);
                 PRECART(A,XS,YS,ZS);
                 IF NOT ECLIPT THEN ECLEQU(TEQX,XS,YS,ZS);
                 IF MODE IN ['G','g']
                   THEN (* -> geocentric *)
                     BEGIN
                       X:=X+XS; Y:=Y+YS; Z:=Z+ZS;
                       WRITELN(' Geocentric coordinates: ');
                     END
                   ELSE (* -> heliocentric *)
                     BEGIN
                       X:=X-XS; Y:=Y-YS; Z:=Z-ZS;
                       WRITELN(' Heliocentric coordinates: ');
                     END;
               END;
       END;

    IF NOT (MODE IN ['?','S','s']) THEN RESULT(X,Y,Z,ECLIPT);

  UNTIL MODE IN ['S','s'];

END. (* COCO *)


