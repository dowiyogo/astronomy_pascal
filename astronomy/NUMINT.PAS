(*-----------------------------------------------------------------------*)
(*                                                                       *)
(*                                NUMINT                                 *)
(*        Numerical integration of perturbed minor planet orbits         *)
(*                           version 93/07/01                            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROGRAM NUMINT ( INPUT,OUTPUT, NUMINP );

  USES  MATLIB,TIMLIB,SPHLIB,PNULIB,KEPLIB,SUNLIB,DELIB,PLALIB;

  CONST J2000  = 0.0;

  VAR DAY,MONTH,YEAR,NLINE,IFLAG         : INTEGER;
      D,HOUR,T_EPOCH, TEQX0,TEQX, FAC    : REAL;
      A,E,M,INC,LAN,AOP                  : REAL;
      MJD_START,MJD_END, T,T1,DT,T2      : REAL;
      XX,YY,ZZ, VX,VY,VZ, XS,YS,ZS       : REAL;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0 : REAL;
      PQR                                : REAL33;
      EQX0_TO_J2000, J2000_TO_EQX, AS    : REAL33;
      Y                                  : DE_EQN_VECTOR;      
      WORK                               : DE_WORKSPACE_RECORD; 
      NUMINP                             : TEXT;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* WRTELM: writes orbital elements                                       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE WRTELM ( YEAR,MONTH: INTEGER; D: REAL;
                   A,E,INC,LAN,AOP,M, TEQX: REAL );
  BEGIN
    WRITELN('  Epoch (y m d)           ',YEAR:7,MONTH:3,D:6:2);
    WRITELN('  Semi-major axis (a)     ',  A:14:7,' AU  ');
    WRITELN('  Eccentricity (e)        ',  E:14:7);
    WRITELN('  Inclination (i)         ',INC:12:5,' deg');
    WRITELN('  Long. of ascending node ',LAN:12:5,' deg');
    WRITELN('  Argument of perihelion  ',AOP:12:5,' deg');
    WRITELN('  Mean anomaly (M)        ',  M:12:5,' deg');
    WRITELN('  Equinox                 ',2000.0+100.0*TEQX:9:2);
  END;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GETELM: reads orbital elements from file NUMINP                       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE GETELM (VAR T,A,E,M: REAL; VAR PQR: REAL33; VAR TEQX0: REAL);

  VAR YEAR,MONTH,DAY : INTEGER;
      D,HOUR         : REAL;
      INC,LAN,AOP    : REAL;

  BEGIN

    WRITELN;
    WRITELN ('                          NUMINT                          ');
    WRITELN ('  Numerical integration of perturbed minor planet orbits  ');
    WRITELN ('                     version 93/07/01                     ');
    WRITELN ('        (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    WRITELN;
    WRITELN (' Orbital elements from file NUMINP: ');
    WRITELN;

    (* Open file for reading *)

    (* RESET(NUMINP); *)                         (* Standard Pascal *)
    ASSIGN(NUMINP,'NUMINP.DAT'); RESET(NUMINP);  (* TURBO Pascal *)

    (* Display orbital elements *)

    READLN (NUMINP,YEAR,MONTH,D);
    READLN (NUMINP,A);
    READLN (NUMINP,E);
    READLN (NUMINP,INC);
    READLN (NUMINP,LAN);
    READLN (NUMINP,AOP);
    READLN (NUMINP,M);
    READLN (NUMINP,TEQX0);

    DAY:=TRUNC(D); HOUR:=24.0*(D-DAY);
    T := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    TEQX0 := (TEQX0-2000.0)/100.0;
    GAUSVEC(LAN,INC,AOP,PQR);

    WRTELM (YEAR,MONTH,D, A,E,INC,LAN,AOP,M, TEQX0 );

  END;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GETEPH: reads desired period of time and equinox of the ephemeris     *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE GETEPH ( VAR T1,DT,T2,TEQX: REAL );

  VAR YEAR,MONTH,DAY :  INTEGER;
      EQX,HOUR,JD :     REAL;

  BEGIN

    WRITELN;
    WRITELN(' Begin and end of the ephemeris: ');
    WRITELN;
    WRITE  ('  First date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  Final date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  Step size (dd hh.hh)                      ... ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
    WRITELN;
    WRITE  (' Desired equinox of the ephemeris (yyyy.y)  ... ');
    READLN (EQX);
    TEQX := (EQX-2000.0)/100.0;

  END;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* WRTLBR: write L and B in deg,min,sec and R                            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE WRTLBR(L,B,R:REAL);
  VAR H,M: INTEGER;
      S  : REAL;
  BEGIN
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); WRITELN(H:5,M:3,TRUNC(S+0.5):3,R:11:6);
  END;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* ACCEL: computes the acceleration                                      *)
(*                                                                       *)
(*   T         Time in Julian centuries since J2000                      *)
(*   X,Y,Z     Heliocentric ecliptic coordinates (in AU)                 *)
(*   AX,AY,AZ  Acceleration (in AU/d**2)                                 *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE ACCEL ( T, X,Y,Z: REAL; VAR AX,AY,AZ: REAL );

  CONST K_GAUSS =  0.01720209895;     (* Gaussian gravitational constant *)

  VAR   PLANET                :  PLANET_TYPE;
        GM_SUN, R_SQR,R,MU_R3 :  REAL;
        XP,YP,ZP, DX,DY,DZ    :  REAL;
        MU                    :  ARRAY[PLANET_TYPE] OF REAL;

  BEGIN

    (* Grav. constant * solar and planetary masses in AU**3/d**2 *)

    GM_SUN := K_GAUSS*K_GAUSS;

    MU[MERCURY] := GM_SUN / 6023600.0;
    MU[VENUS]   := GM_SUN /  408523.5;
    MU[EARTH]   := GM_SUN /  328900.5;
    MU[MARS]    := GM_SUN / 3098710.0;
    MU[JUPITER] := GM_SUN /    1047.355;
    MU[SATURN]  := GM_SUN /    3498.5;
    MU[URANUS]  := GM_SUN /   22869.0;
    MU[NEPTUNE] := GM_SUN /   19314.0;
    MU[PLUTO]   := GM_SUN / 3000000.0;

    (* Solar attraction *)

    R_SQR := ( X*X + Y*Y + Z*Z );  R := SQRT(R_SQR);
    MU_R3 := GM_SUN / (R_SQR*R);
    AX := -MU_R3*X;  AY := -MU_R3*Y;  AZ := -MU_R3*Z;

    (* Planetary perturbation *)

    FOR PLANET := MERCURY TO PLUTO DO

      BEGIN

        (* Planetary coordinates *)

        POSITION ( PLANET, T, XP,YP,ZP );
        DX:=X-XP;  DY:=Y-YP;  DZ:=Z-ZP;

        (* Direct acceleration   *)

        R_SQR := ( DX*DX + DY*DY + DZ*DZ );  R := SQRT(R_SQR);
        MU_R3 := MU[PLANET] / (R_SQR*R);
        AX := AX-MU_R3*DX;  AY := AY-MU_R3*DY;  AZ := AZ-MU_R3*DZ;

        (* Indirect acceleration *)

        R_SQR := ( XP*XP + YP*YP + ZP*ZP );  R := SQRT(R_SQR);
        MU_R3 := MU[PLANET] / (R_SQR*R);
        AX := AX-MU_R3*XP;  AY := AY-MU_R3*YP;  AZ := AZ-MU_R3*ZP;

      END;

  END;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* F: computes the time derivative of the state vector                   *)
(*                                                                       *)
(*   X   Time (Modified Julian Date)                                     *)
(*   Y   State vector (x,y,z in AU, vx,vy,vz in AU/d)                    *)
(*   YP  Derivative (vx,vy,vz in AU/d, ax,ay,az in AU/d**2               *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

{$F+ Turbo Pascal compiler directive: force far call}

PROCEDURE F ( X: REAL; Y: DE_EQN_VECTOR; VAR DYDX: DE_EQN_VECTOR );

  VAR T : REAL;

  BEGIN

    (* Time in Julian centuries since J2000 *)

    T := ( X-51544.5 ) / 36525.0;

    (* Derivative of the state vector *)

    DYDX[1]:=Y[4]; DYDX[2]:=Y[5]; DYDX[3]:=Y[6];

    ACCEL ( T, Y[1],Y[2],Y[3], DYDX[4],DYDX[5],DYDX[6] );

  END;

{$F- Turbo Pascal compiler directive: end far call}


(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* INTEGRATE: Integrates the equation of motion                          *)
(*                                                                       *)
(*   Y        State vector (x,y,z in AU, vx,vy,vz in AU/d)               *)
(*   MJD      Epoch (Modified Julian Date)                               *)
(*   MJD_END  Final epoch (Modified Julian Date)                         *)
(*   IFLAG    Return code                                                *)
(*   WORK     Work space                                                 *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE INTEGRATE ( VAR Y            : DE_EQN_VECTOR;
                      VAR MJD, MJD_END : REAL;
                      VAR IFLAG        : INTEGER;
                      VAR WORK         : DE_WORKSPACE_RECORD );

  CONST EPS = 1.0E-8;  (* Accuracy *)

  VAR   RELERR, ABSERR: REAL;

  BEGIN

    RELERR := EPS;    (* Relative accuracy requirement *)
    ABSERR := 0.0;    (* Absolute accuracy requirement *)

    IF ( MJD_END <> MJD ) THEN
      BEGIN
        REPEAT
          DE ( F, 6, Y, MJD, MJD_END, RELERR,ABSERR, IFLAG, WORK );
        UNTIL ((ABS(IFLAG)=2) OR (IFLAG=6));
        IF (IFLAG=6) THEN WRITELN (' Illegal input in DE ');
      END;

  END;

(*-----------------------------------------------------------------------*)


BEGIN (* NUMINT *)


  (* Read orbital elements, prediction interval and equinox *)

  GETELM (T_EPOCH,A,E,M,PQR,TEQX0);
  GETEPH (T1,DT,T2,TEQX);

  (* Calculate precession matrices *)

  PMATECL (TEQX0,J2000,EQX0_TO_J2000);
  PMATECL (J2000,TEQX,J2000_TO_EQX);

  (* Heliocentric position and velocity vector at epoch *)
  (* referred to ecliptic and equinox of J2000          *)

  ELLIP  ( M,A,E, XX,YY,VX,VY );
  ORBECL ( XX,YY,PQR, Y[1],Y[2],Y[3] );
  ORBECL ( VX,VY,PQR, Y[4],Y[5],Y[6] );

  PRECART ( EQX0_TO_J2000,Y[1],Y[2],Y[3] );
  PRECART ( EQX0_TO_J2000,Y[4],Y[5],Y[6] );


  (* Start integration: propagate state vector from epoch *)
  (* to start of ephemeris                                *)

  MJD_START := T_EPOCH*36525.0 + 51544.5;
  MJD_END   := T1*36525.0 + 51544.5;

  IFLAG     := 1;   (* Initialization flag *)

  INTEGRATE ( Y, MJD_START, MJD_END, IFLAG, WORK );


  (* Orbital elements at start of ephemeris *)

  XX:=Y[1]; YY:=Y[2]; ZZ:=Y[3];    (* Copy J2000 state vector *)
  VX:=Y[4]; VY:=Y[5]; VZ:=Y[6];

  PRECART (J2000_TO_EQX,XX,YY,ZZ);              (* Precession *)
  PRECART (J2000_TO_EQX,VX,VY,VZ);

  XYZKEP ( XX,YY,ZZ,VX,VY,VZ, A,E,INC,LAN,AOP,M ); (* Convert *)

  WRITELN; WRITELN;
  WRITELN (' Orbital elements at start epoch:');
  WRITELN;

  CALDAT (MJD_END,DAY,MONTH,YEAR,HOUR);

  WRTELM ( YEAR,MONTH,DAY+HOUR/24.0,        (* Print elements *)
           A,E,INC,LAN,AOP,M, TEQX );


  (* Create ephemeris *)

  WRITELN; WRITELN;
  WRITELN ('    Date      ET   Sun     l      b     r',
           '        RA          Dec      Distance ');
  WRITELN (' ':45,'   h  m  s      o  ''  "     (AU) ');

  NLINE := 0;
  T     := T1;

  REPEAT

    (* Integrate orbit to time T *)

    MJD_END := T*36525.0 + 51544.5;

    INTEGRATE ( Y, MJD_START, MJD_END, IFLAG, WORK );

    (* Heliocentric ecliptic coordinates, equinox TEQX *)

    XX:=Y[1]; YY:=Y[2]; ZZ:=Y[3];    (* Copy J2000 state vector *)
    VX:=Y[4]; VY:=Y[5]; VZ:=Y[6];

    PRECART (J2000_TO_EQX,XX,YY,ZZ);              (* Precession *)
    PRECART (J2000_TO_EQX,VX,VY,VZ);

    POLAR (XX,YY,ZZ,R,B,L);


    (* Ecliptic coordinates of the Sun, equinox TEXQ *)

    SUN200 (T,LS,BS,RS);  CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS);  PRECART (AS,XS,YS,ZS);


    (* Geometric geocentric coordinates *)

    XX:=XX+XS; YY:=YY+YS; ZZ:=ZZ+ZS;


    (* First order correction for light travel time *)

    DELTA0 := SQRT ( XX*XX + YY*YY + ZZ*ZZ );
    FAC    := 0.00578*DELTA0;

    XX:=XX-FAC*VX;  YY:=YY-FAC*VY;  ZZ:=ZZ-FAC*VZ;

    ECLEQU (TEQX,XX,YY,ZZ);
    POLAR (XX,YY,ZZ,DELTA,DEC,RA); RA:=RA/15.0;


    (* Output *)

    CALDAT (MJD_END,DAY,MONTH,YEAR,HOUR);

    WRITE(YEAR:4,'/',MONTH:2,'/',DAY:2,HOUR:6:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);  WRTLBR(RA,DEC,DELTA0);
    NLINE := NLINE+1; IF (NLINE MOD 5) = 0 THEN WRITELN;


    (* Next time step *)

    T := T + DT;


  UNTIL (T2<T);


END.


