(*-----------------------------------------------------------------------*)
(*                               ORBDET                                  *)
(*         Gaussian orbit determination from three observations          *)
(*               using the abbreviated method of Bucerius                *)
(*                         version 93/07/01                              *)
(*-----------------------------------------------------------------------*)

PROGRAM ORBDET(INPUT,OUTPUT,ORBINP,ORBOUT);

  USES MATLIB, PNULIB, SPHLIB, SUNLIB, TIMLIB, KEPLIB;

  TYPE CHAR80 = ARRAY[1..80] OF CHAR;

  VAR  TEQX                : REAL;
       TP,Q,ECC,INC,LAN,AOP: REAL;
       JD0                 : REAL3;
       RSUN,E              : MAT3X;
       HEADER              : CHAR80;
       ORBINP,ORBOUT       : TEXT;


(*-----------------------------------------------------------------------*)
(* START: reads the input file and preprocesses the observational data   *)
(*                                                                       *)
(* output:                                                               *)
(*   RSUN:  matrix of three Sun position vectors in ecliptic coordinates *)
(*   E:     matrix of three observation direction unit vectors           *)
(*   JD:    julian date of the three observation times                   *)
(*   TEQX:  equinox of RSUN and E (in Julian centuries since J2000)      *)
(*-----------------------------------------------------------------------*)

PROCEDURE START (VAR HEADER: CHAR80;
                 VAR RSUN,E: MAT3X; VAR JD0: REAL3; VAR TEQX: REAL);

  VAR DAY,MONTH,YEAR,D,M,I    :  INTEGER;
      UT,S,DUMMY              :  REAL;
      EQX0,EQX,TEQX0          :  REAL;
      LS,BS,RS,LP,BP,RA,DEC,T :  REAL3;
      A,AS                    :  REAL33;
      ORBINP                  :  TEXT;

  BEGIN

    (* open input file                                                   *)

    (* RESET(ORBINP); *)                              (* Standard Pascal *)
    ASSIGN(ORBINP,'ORBINP.DAT'); RESET(ORBINP);       (* Turbo Pascal    *)

    (* read data from file ORBINP                                        *)

    FOR I:=1 TO 80 DO                                 (* header          *)
      IF NOT(EOLN(ORBINP)) THEN READ(ORBINP,HEADER[I]) ELSE HEADER[I]:=' ';
    READLN(ORBINP);
    FOR I := 1 TO 3 DO                                (* 3 observations  *)
      BEGIN
        READ  (ORBINP,YEAR,MONTH,DAY,UT);                        (* date *)
        READ  (ORBINP,D,M,S); DDD(D,M,S,RA[I]);                  (* RA   *)
        READLN(ORBINP,D,M,S); DDD(D,M,S,DEC[I]);                 (* Dec  *)
        RA[I]:=RA[I]*15.0;
        JD0[I] := 2400000.5+MJD(DAY,MONTH,YEAR,UT);
        T[I]   := (JD0[I]-2451545.0)/36525.0;
      END;
    WRITELN;
    READLN(ORBINP,EQX0); TEQX0:=(EQX0-2000.0)/100.0;  (* equinox         *)

    (* desired equinox of the orbital elements                           *)

    READ(ORBINP,EQX ); TEQX :=(EQX -2000.0)/100.0;

    (* calculate initial data of the orbit determination                 *)

    PMATECL(TEQX0,TEQX,A);
    FOR I := 1 TO 3 DO
      BEGIN
        CART   (1.0,DEC[I],RA[I],E[I,X],E[I,Y],E[I,Z]);
        EQUECL (     TEQX0      ,E[I,X],E[I,Y],E[I,Z]);
        PRECART(      A         ,E[I,X],E[I,Y],E[I,Z]);
        POLAR  (E[I,X],E[I,Y],E[I,Z],DUMMY,BP[I],LP[I]);
        PMATECL( T[I],TEQX,AS);
        SUN200 (T[I],LS[I],BS[I],RS[I]);
        CART   (RS[I],BS[I],LS[I],RSUN[I,X],RSUN[I,Y],RSUN[I,Z]);
        PRECART(      AS         ,RSUN[I,X],RSUN[I,Y],RSUN[I,Z]);
      END;

    WRITELN('   ORBDET: orbit determination from three observations ');
    WRITELN('                     version 93/07/01                  ');
    WRITELN('        (c) 1993 Thomas Pfleger, Oliver Montenbruck    ');
    WRITELN; WRITELN;
    WRITELN(' Summary of orbit determination ');
    WRITELN;
    WRITE  ('  '); FOR I:=1 TO 78 DO WRITE(HEADER[I]); WRITELN;
    WRITELN;
    WRITELN(' Initial data (ecliptic geocentric coordinates (in deg))');
    WRITELN;
    WRITELN('  Julian Date        ', JD0[1]:12:2,JD0[2]:12:2,JD0[3]:12:2);
    WRITELN('  Solar longitude    ',  LS[1]:12:2, LS[2]:12:2, LS[3]:12:2);
    WRITELN('  Planet/Comet Longitude',LP[1]:9:2, LP[2]:12:2, LP[3]:12:2);
    WRITELN('  Planet/Comet Latitude ',BP[1]:9:2, BP[2]:12:2, BP[3]:12:2);
    WRITELN; WRITELN;

  END;

(*-----------------------------------------------------------------------*)
(* DUMPELEM: output of orbital elements (screen)                         *)
(*-----------------------------------------------------------------------*)
PROCEDURE DUMPELEM(TP,Q,ECC,INC,LAN,AOP,TEQX:REAL);
  VAR DAY,MONTH,YEAR: INTEGER;
      MODJD,UT      : REAL;
  BEGIN
    MODJD := TP*36525.0 + 51544.5;
    CALDAT( MODJD, DAY,MONTH,YEAR,UT);
    WRITELN(' Orbital elements',
            ' (Equinox ','J',100.0*TEQX+2000.0:8:2,')');
    WRITELN;
    WRITELN('  Perihelion date      tp    ',
            YEAR:4,'/',MONTH:2,'/',DAY:2,UT:8:4,'h',
            '  (JD',MODJD+2400000.5:11:2,')');
    WRITELN('  Perihelion distance  q[AU] ',  Q:12:6);
    WRITELN('  Semi-major axis      a[AU] ',  Q/(1-ECC):12:6);
    WRITELN('  Eccentricity         e     ',  ECC:12:6);
    WRITELN('  Inclination          i     ',  INC:10:4,' degrees');
    WRITELN('  Ascending node       Omega ',  LAN:10:4,' degrees');
    WRITELN('  Long. of perihelion  pi    ',  AOP+LAN:10:4,' degrees');
    WRITELN('  Arg. of perihelion   omega ',  AOP:10:4,' degrees');
    WRITELN;
  END;

(*-----------------------------------------------------------------------*)
(* SAVEELEM: output of orbital elements (file)                           *)
(*-----------------------------------------------------------------------*)
PROCEDURE SAVEELEM(TP,Q,ECC,INC,LAN,AOP,TEQX:REAL;HEADER: CHAR80);

  VAR I,DAY,MONTH,YEAR: INTEGER;
      MODJD,UT        : REAL;

  BEGIN

    (* open file for writing *)

    (* REWRITE(ORBOUT); *)                            (* Standard Pascal *)
    ASSIGN(ORBOUT,'ORBOUT.DAT'); REWRITE(ORBOUT);     (* Turbo Pascal    *)

    MODJD := TP*36525.0 + 51544.5;
    CALDAT( MODJD, DAY,MONTH,YEAR,UT);
    WRITE  (ORBOUT,YEAR:5,MONTH:3,(DAY+UT/24.0):7:3,'!':6);
    WRITELN(ORBOUT,' perihelion time T0 (y m d.d)  =  JD ',
                   (MODJD+2400000.5):12:3);
    WRITELN(ORBOUT, Q :12:6,'!': 9,' q  ( a =',Q/(1-ECC):10:6,' )');
    WRITELN(ORBOUT,ECC:12:6,'!': 9,' e ');
    WRITELN(ORBOUT,INC:10:4,'!':11,' i ');
    WRITELN(ORBOUT,LAN:10:4,'!':11,' long.asc.node ');
    WRITELN(ORBOUT,AOP:10:4,'!':11,
                   ' arg.perih. ( long.per. = ',AOP+LAN:9:4,' )');
    WRITELN(ORBOUT,TEQX*100.0+2000.0:8:2,'!':13,' equinox (J)');
    WRITE  (ORBOUT,'! ');
    FOR I:=1 TO 78 DO WRITE(ORBOUT,HEADER[I]);

    RESET(ORBOUT); (* close file *)

  END;

(*-----------------------------------------------------------------------*)
(* RETARD: light-time correction                                         *)
(*   JD0: times of observation (t1',t2',t3') (Julian Date)               *)
(*   RHO: three geocentric distances (in AU)                             *)
(*   JD:  times of light emittance (t1,t2,t3) (Julian Date)              *)
(*   TAU: scaled time differences                                        *)
(*-----------------------------------------------------------------------*)
PROCEDURE RETARD ( JD0,RHO: REAL3; VAR JD,TAU: REAL3);
  CONST KGAUSS = 0.01720209895;  A = 0.00578;
  VAR   I: INTEGER;
  BEGIN
    FOR I:=1 TO 3 DO  JD[I]:=JD0[I]-A*RHO[I];
    TAU[1] := KGAUSS*(JD[3]-JD[2]);  TAU[2] := KGAUSS*(JD[3]-JD[1]);
    TAU[3] := KGAUSS*(JD[2]-JD[1]);
  END;

(*-----------------------------------------------------------------------*)
(* GAUSS: iteration of the abbreviated Gauss method                      *)
(*                                                                       *)
(*  RSUN: three vectors of geocentric Sun positions                      *)
(*  E   : three unit vectors of geocentric observation directions        *)
(*  JD0 : three observation times (Julian Date)                          *)
(*  TP  : time of perihelion passage (Julian centuries since J2000)      *)
(*  Q   : perihelion distance                                            *)
(*  ECC : eccentricity                                                   *)
(*  INC : inclination                                                    *)
(*  LAN : longitude of the ascending node                                *)
(*  AOP : argument of perihelion                                         *)
(*-----------------------------------------------------------------------*)

PROCEDURE GAUSS ( RSUN,E: MAT3X; JD0:REAL3;
                  VAR TP,Q,ECC,INC,LAN,AOP: REAL );

  CONST EPS_RHO =1.0E-8;

  VAR I,J              : INTEGER;
      S                : INDEX;
      RHOOLD,DET       : REAL;
      JD,RHO,N,TAU,ETA : REAL3;
      DI               : VECTOR;
      RPL              : MAT3X;
      DD               : REAL33;

  BEGIN

    (* calculate initial approximations of n1 and n3 *)

    N[1] := (JD0[3]-JD0[2]) / (JD0[3]-JD0[1]);     N[2] := -1.0;
    N[3] := (JD0[2]-JD0[1]) / (JD0[3]-JD0[1]);

    (* calculate matrix D and its determinant (det(D) = e3.d3) *)

    CROSS(E[2],E[3],DI);  FOR J:=1 TO 3 DO DD[1,J]:=DOT(DI,RSUN[J]);
    CROSS(E[3],E[1],DI);  FOR J:=1 TO 3 DO DD[2,J]:=DOT(DI,RSUN[J]);
    CROSS(E[1],E[2],DI);  FOR J:=1 TO 3 DO DD[3,J]:=DOT(DI,RSUN[J]);
    DET := DOT(E[3],DI);

    WRITELN; WRITELN(' Iteration of the geocentric distances rho [AU] ');
    WRITELN;

    RHO[2] := 0;

    (* Iterate until distance rho[2] does not change any more *)

    RHO[2] := 0;

    REPEAT

       RHOOLD := RHO[2];

      (* geocentric distance rho from n1 and n3 *)
      FOR I := 1 TO 3 DO
        RHO[I]:=( N[1]*DD[I,1] - DD[I,2] + N[3]*DD[I,3] ) / (N[I]*DET);

      (* apply light-time correction and calculate time differences *)
      RETARD (JD0,RHO,JD,TAU);

      (* heliocentric coordinate vectors *)
      FOR I := 1 TO 3 DO
        FOR S := X TO Z DO
          RPL[I,S] := RHO[I]*E[I,S]-RSUN[I,S];

      (* sector/triangle ratios eta[i] *)
      ETA[1] := FIND_ETA( RPL[2], RPL[3], TAU[1] );
      ETA[2] := FIND_ETA( RPL[1], RPL[3], TAU[2] );
      ETA[3] := FIND_ETA( RPL[1], RPL[2], TAU[3] );

      (* improvement of the sector/triangle ratios *)
      N[1] := ( TAU[1]/ETA[1] ) / (TAU[2]/ETA[2]);
      N[3] := ( TAU[3]/ETA[3] ) / (TAU[2]/ETA[2]);
      WRITELN('  rho',' ':16,RHO[1]:12:8,RHO[2]:12:8,RHO[3]:12:8);

    UNTIL ( ABS(RHO[2]-RHOOLD) < EPS_RHO );

    WRITELN; WRITELN(' Heliocentric distances [AU]:'); WRITELN;
    WRITELN('  r  ',' ':16,
            NORM(RPL[1]):12:8,NORM(RPL[2]):12:8,NORM(RPL[3]):12:8);
    WRITELN; WRITELN;

    (* derive orbital elements from first and third observation *)

    ELEMENT ( JD[1],JD[3],RPL[1],RPL[3], TP,Q,ECC,INC,LAN,AOP );


  END;

(*----------------------------------------------------------------------*)

BEGIN

  START(HEADER,RSUN,E,JD0,TEQX);

  GAUSS(RSUN,E,JD0,TP,Q,ECC,INC,LAN,AOP);

  DUMPELEM(TP,Q,ECC,INC,LAN,AOP,TEQX);
  SAVEELEM(TP,Q,ECC,INC,LAN,AOP,TEQX,HEADER);

  (* check solution  *)

  WRITELN;
  IF (DOT(E[2],RSUN[2])>0) THEN
    WRITELN (' Warning: observation in hemisphere of conjunction;',
             '  possible second solution');
  IF (ECC>1.1) THEN
    WRITELN (' Warning: probably not a realistic solution (e>1.1) ');
  IF ( (ABS(Q-0.985)<0.1) AND (ABS(ECC-0.015)<0.05) ) THEN
    WRITELN (' Warning: probably Earth''s orbit solution');

END.

(*-----------------------------------------------------------------------*)

