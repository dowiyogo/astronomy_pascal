(*---------------------------------------------------------------------------*)
(*                                  PLANRISE                                 *)
(*                planetary and solar rising and setting times               *)
(*                              version 93/07/01                             *)
(*---------------------------------------------------------------------------*)

PROGRAM PLANRISE ( INPUT, OUTPUT );


  USES  MATLIB,TIMLIB,SPHLIB,PNULIB,PLALIB;


  CONST J2000   = 0.0;          (* Standard epoch J2000                  *)
        SID     = 0.9972696;    (* Conversion sidereal/solar time        *)
        SIN_H0P = -9.890038E-3; (* sin(-34'); altitude value for planets *)
        SIN_H0S = -1.45439E-2;  (* sin(-50'); altitude value for the sun *)

  TYPE  EVENT_TYPE = ( RISING, TRANSIT, SETTING );
        STATE_TYPE = ( NEVER_RISES, OK, CIRCUMPOLAR );

  VAR   PLANET                           : PLANET_TYPE;
        EVENT                            : EVENT_TYPE;
        STATE                            : STATE_TYPE;
        MJD0, ZT, ZT0, D_ZT, LST, D_TAU  : REAL;
        LST_0H, SDA, LAMBDA, PHI, ZONE   : REAL;
        T, DEC, RA, SIN_H0               : REAL;
        SIN_PHI, COS_PHI                 : REAL;
        RA_0H, RA_24H, DEC_0H, DEC_24H   : REAL;
        CNT                              : INTEGER;
        PMAT                             : REAL33;


(*---------------------------------------------------------------------------*)
(*                                                                           *)
(* Formatted output                                                          *)
(*                                                                           *)
(*---------------------------------------------------------------------------*)

PROCEDURE WHM(UT:REAL);

  VAR H,M:INTEGER;

  BEGIN
    UT := TRUNC(UT*60.0+0.5)/60.0; (* rounding to 1 min *)
    H:=TRUNC(UT); M:=TRUNC(60.0*(UT-H)+0.5);
    IF H<10 THEN BEGIN WRITE (' 0'); WRITE (H:1,':') END
            ELSE       WRITE (H:3,':');
    IF M<10 THEN BEGIN WRITE ('0'); WRITE (M:1,' ') END
            ELSE       WRITE (M:2,' ')
  END;

(*---------------------------------------------------------------------------*)
(*                                                                           *)
(* Reduce (time) argument X to the interval [0..24]                          *)
(*                                                                           *)
(*---------------------------------------------------------------------------*)

FUNCTION MOD24 (X: REAL): REAL;
  BEGIN
    IF X >= 0.0
      THEN MOD24 := X - 24.0*TRUNC(X/24.0)
      ELSE MOD24 := X - 24.0*TRUNC(X/24.0) + 24.0;
  END;

(*---------------------------------------------------------------------------*)
(*                                                                           *)
(* Read time and geographical coordinates                                    *)
(*                                                                           *)
(*---------------------------------------------------------------------------*)

PROCEDURE GETINPUT (VAR DATE, LAMBDA, PHI, ZONE: REAL);

  VAR D,M,Y: INTEGER;

  BEGIN

    WRITELN;
    WRITELN(' PLANRISE: planetary and solar rising and setting times ');
    WRITELN('                   version 93/07/01                     ');
    WRITELN('      (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    WRITELN;
    WRITE   (' Date (yyyy mm dd)                    ... ');
    READLN (Y,M,D);
    WRITELN;
    WRITE   (' Observing site: longitude (l<0 east) ... '); READLN(LAMBDA);
    WRITE   ('                 latitude             ... '); READLN(PHI);
    WRITE   ('                 local time - UT (h)  ... '); READLN(ZONE);
    WRITELN;

    DATE := MJD(D,M,Y,-ZONE) (* MJD for 0h local time *)

  END; (* GETINPUT *)


(*---------------------------------------------------------------------------*)
(*                                                                           *)
(* PLAN_RA_DEC: geocentric equatorial planetary coordinates                  *)
(*                                                                           *)
(* PLANET: Planet for which the coordinates are computed. A call with        *)
(*         PLANET = EARTH yields the geocentric coordinartes of the sun.     *)
(* T     : time in julian centuries since J2000                              *)
(* RA    : Right ascension in deg [0..360]                                   *)
(* DEC   : Declination in deg                                                *)
(*                                                                           *)
(* Note:                                                                     *)
(* This procedure uses the globally defined precesion matrix PMAT            *)
(* for transformation from J2000 to the required equinox.                    *)
(*                                                                           *)
(*---------------------------------------------------------------------------*)

PROCEDURE PLAN_RA_DEC (PLANET: PLANET_TYPE; T: REAL; VAR RA, DEC: REAL);

  VAR XP,YP,ZP, XE,YE,ZE, R: REAL;

  BEGIN

  IF (PLANET<>EARTH) THEN

    BEGIN

      POSITION ( EARTH, T, XE,YE,ZE );
      ECLEQU   ( J2000,    XE,YE,ZE );

      (* Determine geocentric geometric planetary position *)

      POSITION ( PLANET, T, XP,YP,ZP );
      ECLEQU   ( J2000, XP,YP,ZP );

      XP := XP-XE;  YP := YP-YE;  ZP := ZP-ZE;
 
      PRECART  ( PMAT,  XP,YP,ZP);

      (* Right ascension, declination and distance of the planet *)
      POLAR ( XP,YP,ZP, R, DEC, RA );

    END

  ELSE

    (* Compute geocentric equatorial coordinates of the sun *)

    BEGIN

      POSITION ( EARTH, T, XE,YE,ZE );
      ECLEQU   ( J2000, XE,YE,ZE );
      PRECART  ( PMAT,  XE,YE,ZE);

      POLAR ( -XE,-YE,-ZE, R, DEC, RA );

    END;

END; (* PLAN_RA_DEC *)

(*---------------------------------------------------------------------------*)


BEGIN (* PLANRISE main program *)

  (* Read input data *)

  GETINPUT (MJD0, LAMBDA, PHI, ZONE);

  SIN_PHI := SN (PHI);
  COS_PHI := CS (PHI);

  (* Compute local sidereal time at 0h local time *)

  LST_0H := LMST (MJD0,LAMBDA);

  (* Precession matrix (J2000 -> mean equator and equinox of date)  *)
  (* for equatorial coordinates                                     *)

  T := (MJD0 - 51544.5)/36525.0;
  PMATEQU (J2000, T, PMAT);

  (* Compute and print rising and setting times *)

  WRITELN ('rise':18,'culmination':16,'set':9); WRITELN;

  FOR PLANET:=MERCURY TO PLUTO DO

    BEGIN

      (* Compute geocentr. planetary position at 0h and 24h local time *)
      T := (MJD0 - 51544.5)/36525.0;
      PLAN_RA_DEC (PLANET,T,  RA_0H,  DEC_0H);
      T := (MJD0 + 1.0 - 51544.5)/36525.0;
      PLAN_RA_DEC (PLANET,T, RA_24H, DEC_24H);

      (* Generate continuous right ascension values in case of jumps *)
      (* between 0h and 24h                                          *)
      IF (RA_0H - RA_24H) > 180.0 THEN RA_24H := RA_24H + 360.0;
      IF (RA_0H - RA_24H) <-180.0 THEN RA_0H  := RA_0H  + 360.0;

      CASE PLANET OF
        MERCURY: WRITE (' Mercury ');
        VENUS  : WRITE (' Venus   ');
        EARTH  : WRITE (' Sun     ');
        MARS   : WRITE (' Mars    ');
        JUPITER: WRITE (' Jupiter ');
        SATURN : WRITE (' Saturn  ');
        URANUS : WRITE (' Uranus  ');
        NEPTUNE: WRITE (' Neptune ');
        PLUTO  : WRITE (' Pluto   ')
        END;
      WRITE(' ':3);

      EVENT := RISING; STATE := OK;

      WHILE ( (EVENT<=SETTING) AND (STATE=OK) ) DO

        BEGIN

          ZT0 := 12.0; (* Starting value 12h local time *)
          CNT := 0;

          REPEAT

            (* Linear interpolation of planetary position *)
            
            RA  := RA_0H  + (ZT0/24.0) * (RA_24H - RA_0H);
            DEC := DEC_0H + (ZT0/24.0) * (DEC_24H - DEC_0H);

            (* Compute semi-diurnal arc (in deg) *)

            IF PLANET<>EARTH
              THEN SIN_H0 := SIN_H0P
              ELSE SIN_H0 := SIN_H0S;

            SDA := (SIN_H0 - SN (DEC) * SIN_PHI ) / (CS (DEC) * COS_PHI);

            IF (ABS(SDA)<1.0)
              THEN 
                BEGIN 
                  SDA  :=ACS(SDA);
                  STATE:=OK
                END
              ELSE (* Test for circumpolar motion or invisibility *)
                IF  (PHI>=0.0) 
                  THEN
                    IF DEC>(90.0-PHI)
                      THEN STATE:=CIRCUMPOLAR ELSE STATE:=NEVER_RISES
                  ELSE
                    IF DEC<(-90.0-PHI)
                      THEN STATE:=CIRCUMPOLAR ELSE STATE:=NEVER_RISES;

            (* Improved times for rising, culmination and setting *)

            IF STATE=OK THEN
              BEGIN
                LST := LST_0H+ZT0/SID; (* Sidereal time at univ. time ZT0 *)
                CASE EVENT OF
                  RISING : D_TAU := (LST-RA/15.0) + SDA/15.0;
                  TRANSIT: D_TAU := (LST-RA/15.0);
                  SETTING: D_TAU := (LST-RA/15.0) - SDA/15.0;
                END;
                D_ZT := SID*(MOD24(D_TAU+12.0)-12.0);
                ZT := ZT0-D_ZT;
                ZT0 := ZT;
                CNT := CNT+1
              END;


          UNTIL ( (ABS(D_ZT)<=0.008) OR (CNT>10) OR (STATE<>OK));


          (* Print result *)

          IF STATE=OK
          THEN BEGIN WHM (ZT); WRITE (' ':6); END
          ELSE
            CASE STATE OF
              NEVER_RISES: WRITE ('-------- always invisible -------');
              CIRCUMPOLAR: WRITE ('--------- always visible ---------')
            END;

          EVENT := SUCC (EVENT);

      END;

      WRITELN;

    END; (* FOR PLANET... *)

  WRITELN;  WRITE  (' all times in local time ( = UT ');
  IF ZONE>=0 THEN WRITE('+'); WRITELN(ZONE:4:1,'h )');

END.

