(*-----------------------------------------------------------------------*)
(*                                PHYS                                   *)
(*        physical ephemerides of the major planets and the Sun          *)
(*                         version 93/07/01                              *)
(*-----------------------------------------------------------------------*)

PROGRAM PHYS ( INPUT, OUTPUT );

  USES  MATLIB,TIMLIB,SPHLIB,PNULIB,PLALIB,PHYLIB;

  CONST AU      = 149597870.0; (* 1 AU in km                 *)
        C_LIGHT = 173.14;      (* Velocity of light ( AU/d ) *)
        J2000   = 0.0;         (* Reference epoch J2000      *)

  VAR   PLANET               : PLANET_TYPE;
        DAY,MONTH,YEAR       : INTEGER;
        HOUR, T,T0           : REAL;
        X,Y,Z, XE,YE,ZE      : REAL;
        XX,YY,ZZ, R, DELTA   : REAL;
        R_EQU, F, D_EQU      : REAL;
        DD, DM               : INTEGER;
        DS                   : REAL;
        L1,L2,L3, B, BSUN, D : REAL;
        A0,D0,W1,W2,W3       : REAL;
        SENSE                : ROTATION_TYPE;
        AX,AY,AZ, POSAX      : REAL;
        ELONG, PHI,K, MAG    : REAL;
        LSUN, DSUN, POSSUN   : REAL;
        PMAT                 : REAL33;


BEGIN

   (* Print header and read desired date *)

  WRITELN;
  WRITELN (' PHYS: physical ephemerides of the planets and the Sun');
  WRITELN ('                 Version 93/07/01                     ');
  WRITELN ('     (c) 1993 Thomas Pfleger, Oliver Montenbruck      ');
  WRITELN;
  WRITE   (' Date (yyyy mm dd hh.hhh)  ');
  READLN  (YEAR,MONTH,DAY,HOUR);
  WRITELN;

  T := ( MJD(DAY,MONTH,YEAR,HOUR)-51544.5 ) / 36525.0;

  WRITE   ('D':13,'V':7,'i':7,'PA(S)':9,'PA(A)':8);
  WRITELN ('L(I)':8,'L(II)':7,'L(III)':7,'B':6);
  WRITELN ('"':13,'mag':8,'o':6,'o':7,'o':8,'o':8,'o':7,'o':7,'o':8);

  (* Precession matrix (J2000 -> mean equinox of date) *)
  (* for equatorial coordinates                        *)

  PMATEQU (J2000, T, PMAT);

  (* Equatorial coordinates of the Earth, mean equinox of J2000 *)

  POSITION ( EARTH, T, XE,YE,ZE );
  ECLEQU   ( J2000, XE,YE,ZE );


  (* Physical ephemerides of the planets *)

  FOR PLANET:=MERCURY TO PLUTO DO
   IF PLANET<>EARTH THEN

    BEGIN

      (* Compute the planet's geocentric geometric position and the *)
      (* light time (in days)                                       *)

      POSITION ( PLANET, T, X,Y,Z );
      ECLEQU   ( J2000, X,Y,Z );
      DELTA := SQRT ( (X-XE)*(X-XE) + (Y-YE)*(Y-YE) + (Z-ZE)*(Z-ZE) );
      T0 := T - (DELTA/C_LIGHT) / 36525.0;

      (* Compute the antedated planetary position at emission of light *)
      (* (i.e. apply a first order light time correction)              *)

      POSITION ( PLANET, T0, X,Y,Z );
      ECLEQU   ( J2000, X,Y,Z );

      (* Light time corrected geocentric coordinates *)

      XX := X-XE;  YY := Y-YE;  ZZ := Z-ZE;

      (* Compute apparent equatorial diameter (in ") *)

      SHAPE ( PLANET, R_EQU, F );
      D_EQU := 3600.0 * 2.0*ASN(R_EQU/(DELTA*AU));

      (* Compute right ascension and declination of the axis of     *)
      (* rotation with respect to the equator and equinox of J2000; *)
      (* compute orientation of the prime meridian                  *)

      ORIENT ( PLANET,SYS_I  ,T0, A0,D0,W1,SENSE );
      ORIENT ( PLANET,SYS_II ,T0, A0,D0,W2,SENSE );
      ORIENT ( PLANET,SYS_III,T0, A0,D0,W3,SENSE );

      (* Compute planetocentric longitude and latitude of the Earth *)

      ROTATION ( XX,YY,ZZ,A0,D0,W1,SENSE,F, AX,AY,AZ,L1,B,D );
      ROTATION ( XX,YY,ZZ,A0,D0,W2,SENSE,F, AX,AY,AZ,L2,B,D );
      ROTATION ( XX,YY,ZZ,A0,D0,W3,SENSE,F, AX,AY,AZ,L3,B,D );

      (* Compute planetocentric longitude and latitude of the Sun *)

      ROTATION ( X,Y,Z,A0,D0,W1,SENSE,F, AX,AY,AZ,LSUN,BSUN,DSUN );

      (* Compute illumination and apparent magnitude *)

      ILLUM ( X,Y,Z, XE,YE,ZE, R,DELTA,ELONG, PHI,K );
      MAG := BRIGHT ( PLANET, R,DELTA,PHI,DSUN,LSUN-L1 );

      (* Compute position angles of the axis of rotation and of the *)
      (* Sun with respect to the mean equinox of date               *)

      PRECART ( PMAT,  X, Y, Z );
      PRECART ( PMAT, XX,YY,ZZ );
      PRECART ( PMAT, AX,AY,AZ );

      POSAX  := POSANG ( XX,YY,ZZ, AX,AY,AZ );
      POSSUN := POSANG ( XX,YY,ZZ, -X,-Y,-Z );

      (* Print results *)

      CASE PLANET OF
        MERCURY: WRITE (' Mercury  ');
        VENUS  : WRITE (' Venus    ');
        MARS   : WRITE (' Mars     ');
        JUPITER: WRITE (' Jupiter  ');
        SATURN : WRITE (' Saturn   ');
        URANUS : WRITE (' Uranus   ');
        NEPTUNE: WRITE (' Neptune  ');
        PLUTO  : WRITE (' Pluto    ')
      END;

      WRITE ( D_EQU:5:2, MAG:6:1, PHI:7:1 );
      WRITE ( POSSUN:8:2, POSAX:8:2 );

      CASE PLANET OF
        MERCURY,VENUS,MARS,PLUTO:  WRITE (L1:8:2, ' ':14);
        JUPITER:                   WRITE (L1:8:2, L2:7:2, L3:7:2);
        SATURN:                    WRITE (L1:8:2, L3:14:2);
        URANUS,NEPTUNE:            WRITE (L3:22:2)
      END;

      WRITELN (B:8:2);

    END;


  (* Physical ephemerides of the Sun  *)

  (* Compute light time corrected equatorial coordinates of     *)
  (* the Earth with respect to the equator and equinox of J2000 *)

  DELTA := SQRT ( XE*XE + YE*YE + ZE*ZE );   
  T0 := T - (DELTA/C_LIGHT) / 36525.0;       
  POSITION ( EARTH, T0, XE,YE,ZE );
  ECLEQU   ( J2000, XE,YE,ZE );

  (* Right ascension and declination of the Sun's axis (J2000),    *)
  (* orientation of the prime meridian and equatorial radius (km)  *)

  A0 := 285.96;  
  D0 :=  63.96;
  W1 := 84.11 + 14.1844000 * 36525.0*T;
  W1 := W1/360.0; W1:=360.0*(W1-TRUNC(W1));
  R_EQU := 696000.0;

  (* Compute heliographic coordinates of the Earth *)

  ROTATION ( -XE,-YE,-ZE, A0,D0,W1, RETROGRADE, 0.0, AX,AY,AZ,L1,B,B );

  (* Compute position angle of the axis of rotation *)
  (* with respect to the mean equinox of date       *)

  PRECART ( PMAT, XE,YE,ZE );
  PRECART ( PMAT, AX,AY,AZ );
  POSAX  := POSANG ( -XE,-YE,-ZE, AX,AY,AZ );

  (* Express position angle of the Sun's axis within -180..180 deg *)

  IF POSAX>180.0 THEN POSAX:=POSAX-360.0;

  (* Compute apparent equatorial diameter (in ") *)

  D_EQU := 2.0*ASN(R_EQU/(DELTA*AU));  

  (* Print results *)

  DMS (D_EQU, DD,DM,DS);
  WRITELN ('''':10,'"':3,'o':29,'o':8,'o':22);
  WRITELN (' Sun  ', DM:3, DS:6:2, POSAX:29:2, L1:8:2, B:22:2 )

END.

