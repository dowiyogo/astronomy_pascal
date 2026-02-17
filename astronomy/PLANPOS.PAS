(*-----------------------------------------------------------------------*)
(*                              PLANPOS                                  *)
(*           heliocentric and geocentric planetary positions             *)
(*                         version 93/07/01                              *)
(*-----------------------------------------------------------------------*)

PROGRAM PLANPOS(INPUT,OUTPUT);

  uses matlib, pnulib, sphlib, sunlib, p15lib, p69lib, timlib;

  CONST J2000 =  0.0;
        B1950 = -0.500002108;

  VAR DAY, MONTH, YEAR, IPLAN, IMODE, K    : INTEGER;
      HOUR, MODJD, T, TEQX                 : REAL;
      X,Y,Z, XP,YP,ZP, XS,YS,ZS            : REAL;
      L,B,R, LS,BS,RS, RA,DEC,DELTA,DELTA0 : REAL;
      A                                    : REAL33;
      MODE                                 : CHAR;

(*-----------------------------------------------------------------------*)
(* PRINTOUT: print coordinates of one planet                             *)
(*-----------------------------------------------------------------------*)
PROCEDURE PRINTOUT(IPLAN:INTEGER;L,B,R,RA,DEC,DELTA:REAL);
  VAR H,M: INTEGER;
      S  : REAL;
  BEGIN
    DMS(L,H,M,S);   WRITE (H:3,M:3,S:5:1);
    DMS(B,H,M,S);   WRITE (H:4,M:3,S:5:1);
    IF IPLAN<4 THEN WRITE (R:11:6)
               ELSE WRITE (R:10:5,' ');
    DMS(RA,H,M,S);  WRITE (H:4,M:3,S:6:2);
    DMS(DEC,H,M,S); WRITE (H:4,M:3,S:5:1);
    IF IPLAN<4 THEN WRITE (DELTA:11:6)
               ELSE WRITE (DELTA:10:5,' ');
  END;
(*-----------------------------------------------------------------------*)

BEGIN (* PLANPOS *)

  WRITELN;
  WRITELN('   PLANPOS: geocentric and heliocentric planetary positions ');
  WRITELN('                       version 93/07/01                     ');
  WRITELN('         (c) 1993 Thomas Pfleger, Oliver Montenbruck        ');
  WRITELN;

  REPEAT

    WRITELN;
    WRITELN (' (J) J2000 astrometric       (B) B1950 astrometric  ');
    WRITELN (' (A) apparent coordinates    (E) end                ');
    WRITELN;
    WRITE   (' enter option: '); READLN (MODE);
    WRITELN;

    IF MODE IN ['A','a','J','j','B','b']  THEN

      BEGIN

        (* read date *)

        WRITE (' date (year month day hour) ?     ');
        READLN (YEAR,MONTH,DAY,HOUR);  WRITELN; WRITELN; WRITELN; WRITELN;
        MODJD := MJD (DAY,MONTH,YEAR,HOUR);  T:=(MODJD-51544.5)/36525.0;
        WRITE (' date:  ', YEAR:4,'/',MONTH:2,'/',DAY:2,' ',HOUR:5:1,'(ET)');
        WRITE ('JD:':6,(MODJD+2400000.5):12:3,'equinox ':18);
        CASE MODE OF
          'A','a': WRITELN ('of date');
          'J','j': WRITELN ('J2000');
          'B','b': WRITELN ('B1950');
         END;
        WRITELN;

        (* header *)

        WRITE   (' ':10,'l':6,'b':12,'r':11);
        WRITELN (' ':7,'RA':5,'Dec':13,'delta':13);
        WRITE   (' ':9,'   o  ''  "',' ':3,'  o  ''  "',' ':6,'AU',' ':4);
        WRITELN (' ':2,'  h  m  s',' ':4,'  o  ''  "',' ':6,'AU');

        (* ecliptic coordinates of the sun, Equinox T  *)

        SUN200 (T,LS,BS,RS);

        (* planetary coordinates; include Pluto between 1890 and 2100  *)

        IF ( (-1.1<T) AND (T<+1.0) ) THEN K:=9 ELSE K:=8;

        FOR IPLAN:=0 TO K DO

          BEGIN

            (* heliocentric ecliptic coordinates of the planet        *)

            CASE IPLAN OF
              1: MER200(T,L,B,R); 2: VEN200(T,L,B,R);
              4: MAR200(T,L,B,R); 5: JUP200(T,L,B,R); 6: SAT200(T,L,B,R);
              7: URA200(T,L,B,R); 8: NEP200(T,L,B,R); 9: PLU200(T,L,B,R);
              0: BEGIN L:=0.0; B:=0.0; R:=0.0; END;
              3: BEGIN L:=LS+180.0; B:=-BS; R:=RS; END;
             END;

            (* geocentric ecliptic coordinates (light-time corrected)  *)

            IF MODE IN ['A','a'] THEN IMODE:=2 ELSE IMODE:=1;
            GEOCEN (T, L,B,R, LS,BS,RS, IPLAN,IMODE,
                    XP,YP,ZP, XS,YS,ZS, X,Y,Z,DELTA0);

            (* precession, equatorial coordinates, nutation            *)

            CASE MODE OF
              'J','j': TEQX:=J2000; 'B','b': TEQX:=B1950;
             END;

            IF MODE IN ['A','a']
              THEN
                BEGIN  ECLEQU(T,X,Y,Z); NUTEQU(T,X,Y,Z); END
              ELSE
                BEGIN
                  PMATECL(T,TEQX,A); PRECART (A,XP,YP,ZP);
                  PRECART (A,X,Y,Z); ECLEQU (TEQX,X,Y,Z);
                END;

            (* spherical coordinates *)

            POLAR (XP,YP,ZP,R,B,L);
            POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

            (* output *)

            CASE IPLAN OF
              0: WRITE(' Sun     ');  1: WRITE(' Mercury ');
              2: WRITE(' Venus   ');  3: WRITE(' Earth   ');
              4: WRITE(' Mars    ');  5: WRITE(' Jupiter ');
              6: WRITE(' Saturn  ');  7: WRITE(' Uranus  ');
              8: WRITE(' Neptune ');  9: WRITE(' Pluto   ');
             END;

            PRINTOUT(IPLAN,L,B,R,RA,DEC,DELTA0); WRITELN;

          END;

          WRITELN;
          WRITELN (' l,b,r:   heliocentric ecliptic (geometric) ');
          WRITE   (' RA,Dec:  geocentric equatorial ');
          IF MODE IN ['A','a'] THEN WRITELN('(apparent)')
                               ELSE WRITELN('(astrometric)');
          WRITELN (' delta:   geocentric distance   (geometric)');
          WRITELN;

      END;

  UNTIL MODE IN ['E','e']

END. (* PLANPOS *)


