(*-----------------------------------------------------------------------*)
(*                                OCCULT                                 *)
(*            prediction of stellar occultations by the Moon             *)
(*                           version 93/07/01                            *)
(*-----------------------------------------------------------------------*)

PROGRAM OCCULT(INPUT,OUTPUT,OCCINP);

  USES MATLIB, PNULIB, SPHLIB, SUNLIB, MOOLIB, TIMLIB;

  CONST TOVLAP      = 3.42E-6;          (*  3h in julian centuries       *)
        T_SEARCH    = 2.737850787E-4;   (* 10d in julian centuries       *)
        NAME_LENGTH = 17;               (* maxim.length of a star's name *)

  TYPE  NAME_STRING  = ARRAY[1..NAME_LENGTH] OF CHAR;

  VAR   T_BEGIN,T_END,T1,T2,TM,T_EQX,T_EPOCH : REAL;
        ETDIFUT                              : REAL;
        RA_STAR,DE_STAR                      : REAL;
        VX,VY,VZ                             : REAL;
        LAMBDA,PHI,RCPHI,RSPHI               : REAL;
        RAPOLY,DEPOLY,RPOLY                  : TPOLYNOM;
        PNMAT                                : REAL33;
        OCCINP                               : TEXT;
        NAME                                 : NAME_STRING;


(*-----------------------------------------------------------------------*)
(* GET_INPUT: read desired period of time and observer's coordinates     *)
(*-----------------------------------------------------------------------*)
PROCEDURE GET_INPUT(VAR T_BEGIN,T_END,ETDIFUT,LAMBDA,PHI: REAL);
  VAR D,M,Y: INTEGER;
      T    : REAL;
      VALID: BOOLEAN;
  BEGIN
    WRITELN;
    WRITELN ('          OCCULT: occultations of stars by the Moon       ');
    WRITELN ('                      version 93/07/01                    ');
    WRITELN ('        (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    WRITELN;
    WRITELN (' Period of time for prediction of occultations ');
    WRITE ('   first date (yyyy mm dd)                   ... ');
    READLN(Y,M,D);
    T_BEGIN := (MJD(D,M,Y,0)-51544.5)/36525.0;
    WRITE ('   last  date (yyyy mm dd)                   ... ');
    READLN(Y,M,D);
    T_END   := (MJD(D,M,Y,0)-51544.5)/36525.0;
    T := ( T_BEGIN + T_END ) / 2.0;
    ETMINUT ( T, ETDIFUT, VALID);
    IF (VALID)
      THEN WRITE(' Difference ET-UT (proposal:',
                  TRUNC(ETDIFUT+0.5):3,' sec)          ... ')
      ELSE WRITE(' Difference ET-UT (sec)                      ... ');
    READLN(ETDIFUT);
    WRITE (' Observer''s coordinates: longitude (>0 west) ... ');
    READLN(LAMBDA);
    WRITE ('                         latitude            ... ');
    READLN(PHI);
  END;

(*-----------------------------------------------------------------------*)
(* HEADER: print header                                                  *)
(*-----------------------------------------------------------------------*)
PROCEDURE HEADER;
  BEGIN
    WRITELN;
    WRITELN ('   Date        UT     D/R   Pos ',
             '    h      a     b       Star');
    WRITELN ('              h  m  s          o',
             '     o     m     m           ');
  END;

(*-----------------------------------------------------------------------*)
(* GETSTAR: read star coordinates from file OCCINP and correct for       *)
(*          proper motion                                                *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETSTAR ( T_EPOCH,T: REAL;
                    VAR RA,DEC: REAL; VAR NAME: NAME_STRING);
  VAR G,M,I            :  INTEGER;
      S, PM_RA, PM_DEC :  REAL;
  BEGIN
    READ(OCCINP,G,M,S);  DDD(G,M,S,RA);      (* right ascension at epoch *)
    READ(OCCINP,PM_RA );
    READ(OCCINP,G,M,S);  DDD(G,M,S,DEC);     (* declination at epoch     *)
    READ(OCCINP,PM_DEC);
    RA  := RA  + (T-T_EPOCH)*PM_RA /3600.0;  (* proper motion right asc. *)
    DEC := DEC + (T-T_EPOCH)*PM_DEC/3600.0;  (* proper motion declination*)
    FOR I:=1 TO NAME_LENGTH DO               (* name of the star         *)
      IF (NOT EOLN(OCCINP) ) THEN READ(OCCINP,NAME[I]) ELSE NAME[I]:=' ';
    READLN(OCCINP);
    RA := 15.0 * RA;                         (* RA in deg                *)
  END;

(*-----------------------------------------------------------------------*)
(*  CONJUNCT:                                                            *)
(*                                                                       *)
(*  checks whether there is a conjunction of Moon and star between TA    *)
(*  and TB during which the Moon's shadow hits the Earth                 *)
(*                                                                       *)
(*  TA,TB:    time interval for search of conjunction                    *)
(*  RAPOLY,DEPOLY,RPOLY: Chebyshev coefficients for lunar coordinates    *)
(*  RA,DEC:   right ascension and declination of the star (0<=RA<=360)   *)
(*  CONJ:     TRUE/FALSE (conjunction found / no conjunction found)      *)
(*  T_CONJ:   time of conjunction in right ascension (0.0 if CONJ=FALSE) *)
(*                                                                       *)
(*  All times are counted in Julian centuries since J2000.               *)
(*  The Chebyshev expansion of the lunar right ascension has to yield    *)
(*  values between -360 and +360 degress and cover less than one orbit.  *)
(*-----------------------------------------------------------------------*)

PROCEDURE CONJUNCT ( TA,TB: REAL; RAPOLY,DEPOLY,RPOLY: TPOLYNOM;
                     RA,DEC: REAL; VAR CONJ: BOOLEAN; VAR T_CONJ: REAL);

  CONST EPS=1E-3; (* accuracy in degrees RA *)

  VAR RA_A,RA_B                    : REAL;
      T1,T2,T_NEW,DRA1,DRA2,DRA_NEW: REAL;
      DE_CON, R_CON                : REAL;

  BEGIN

    T_CONJ := 0.0;
    RA_A := T_EVAL(RAPOLY,TA);
    RA_B := T_EVAL(RAPOLY,TB);

    (* check if RA_A <= RA <= RA_B *)
    CONJ := (RA_A<=RA) AND (RA<=RA_B);
    IF (NOT CONJ) THEN   (* check again with RA-360deg *)
      BEGIN RA:=RA-360.0; CONJ := ((RA_A<=RA) AND (RA<=RA_B)); END;

    IF CONJ THEN

      BEGIN

        (* determine time of conjunction using 'regula falsi' *)
        (* ([T1,T2] always contains T_CONJ)                   *)

        T1 := TA; DRA1 :=RA_A-RA;
        T2 := TB; DRA2 :=RA_B-RA;
        REPEAT
          T_NEW   := T2 - DRA2*(T2-T1)/(DRA2-DRA1);
          DRA_NEW := T_EVAL(RAPOLY,T_NEW) - RA;
          IF DRA_NEW>0 THEN  BEGIN T2:=T_NEW; DRA2:=DRA_NEW END
                       ELSE  BEGIN T1:=T_NEW; DRA1:=DRA_NEW END;
        UNTIL (ABS(DRA_NEW)<EPS);
        T_CONJ := T_NEW;

        (* check if lunar shadow hits the Earth               *)

        DE_CON := T_EVAL(DEPOLY,T_CONJ);
        R_CON  := T_EVAL(RPOLY, T_CONJ);
        CONJ :=  ( ABS(SN(DE_CON-DEC)*R_CON) < 1.5 );

      END;

  END;

(*-----------------------------------------------------------------------*)
(*  SHADOW:                                                              *)
(*                                                                       *)
(*  starting from the time of conjunction the times, position angles     *)
(*  and longitude and latitude coefficients of disappearance and         *)
(*  reappearance are calculated for a specific observing site            *)
(*                                                                       *)
(*  RAPOLY,DEPOLY,RPOLY  : Chebyshev approximations of lunar coordinates *)
(*  T_CONJ_ET            : time of conjunction in right ascension        *)
(*                         (in Julian centuries ET since J2000)          *)
(*  ETDIFUT              : ET-UT in sec                                  *)
(*  LAMBDA,RCPHI,RSPHI   : geocentric coordinates of the observer        *)
(*  RA_STAR,DE_STAR      : right ascension and declination of the star   *)
(*  EVENT                : TRUE occultation takes place; FALSE otherwise *)
(*  MJD_UT_IN, MJD_UT_OUT: times of contact (Modified Julian Date UT)    *)
(*  POS_IN, POS_OUT      : position angles                               *)
(*  H_IN, H_OUT          : star's altitude above the horizon             *)
(*  A_IN, A_OUT          : longitude coefficient                         *)
(*  B_IN, B_OUT          : latitude coefficient                          *)
(*-----------------------------------------------------------------------*)

PROCEDURE SHADOW

  ( RAPOLY, DEPOLY, RPOLY                                  : TPOLYNOM;
    T_CONJ_ET,ETDIFUT, LAMBDA,RCPHI,RSPHI,  RA_STAR,DE_STAR: REAL;
    VAR EVENT                                              : BOOLEAN;
    VAR MJD_UT_IN, MJD_UT_OUT                              : REAL;
    VAR POS_IN,POS_OUT, H_IN,H_OUT, A_IN,A_OUT, B_IN,B_OUT : REAL      );


  CONST DTAB     = 0.25;   (* search step size in hours                  *)
        RANGE    = 2.25;   (* search interval = +/-(RANGE+DTAB) in hours *)
        K        = 0.2725;          (* ratio earth radius / lunar radius *)
        CENT     = 876600.0;        (* hours per Julian century          *)
        SID      = 1.0027379;       (* ratio solar time / sidereal time  *)

  VAR I, NZ, NFOUND                  : INTEGER;
      K_SQR, MJD_CONJ_UT, HOUR, F, G : REAL;
      THETA_CONJ, CSDEST, SNDEST     : REAL;
      S_MINUS, S_0, S_PLUS, XE, YE   : REAL;
      Z, TIME                        : ARRAY[1..2] OF REAL;

  (* FG: f-g coordinates in the fundamental plane *)
  PROCEDURE FG ( HOUR: REAL; VAR F,G: REAL);
    VAR T,DEM,RM,RCDEM,RSDEM,DRAM,DRA: REAL;
    BEGIN
      T := T_CONJ_ET + HOUR/CENT;
      DEM   := T_EVAL(DEPOLY,T);   RM    := T_EVAL(RPOLY,T);
      RCDEM := RM * CS(DEM);       RSDEM := RM * SN(DEM);
      DRAM  := T_EVAL(RAPOLY,T) - RA_STAR;
      DRA   := 15.0 * ( THETA_CONJ + HOUR*SID ) - RA_STAR;
      F := +RCDEM*SN(DRAM) - RCPHI*SN(DRA) ;
      G := + RSDEM*CSDEST - RCDEM*SNDEST*CS(DRAM)
           - RSPHI*CSDEST + RCPHI*SNDEST*CS(DRA);
    END;

  (* CONTACT: position angle, altitude and longitude/latitude coeffic.   *)
  PROCEDURE CONTACT ( HOUR: REAL; VAR POS,H,A,B: REAL );
    VAR F,G,FF,GG,DF,DG,FAC,DRA,CDRA,SDRA: REAL;
    BEGIN
      FG ( HOUR, F, G );  FG ( HOUR+DTAB , FF, GG );
      DF := (FF-F)/DTAB;  DG := (GG-G)/DTAB;
      POS := ATN2(-F,-G);  IF POS<0.0 THEN POS:=POS+360.0;
      FAC := 1.047 / (F*DF+G*DG);
      DRA := 15.0 * ( THETA_CONJ + HOUR*SID ) - RA_STAR;
      CDRA := CS(DRA); SDRA := SN(DRA);
      A := -FAC * RCPHI * ( F*CDRA + G*SDRA*SNDEST );
      B := -FAC*( RSPHI * (F*SDRA-G*SNDEST*CDRA) - RCPHI*G*CSDEST );
      H :=  ASN ( RSPHI*SNDEST + RCPHI*CSDEST*CDRA );
    END;


  BEGIN

    (* modified julian date and sidereal time at time of conjunction     *)
    MJD_CONJ_UT  := T_CONJ_ET*36525.0 + 51544.5 - ETDIFUT/86400.0;
    THETA_CONJ   := LMST ( MJD_CONJ_UT, LAMBDA );

    (* auxiliary values *)
    K_SQR := K*K;  CSDEST := CS(DE_STAR); SNDEST := SN(DE_STAR);

    (* search for time of contact *)
    NFOUND := 0;  TIME[1]:=0.0; TIME[2]:=0.0;
    HOUR := -RANGE-2.0*DTAB;
    FG (-RANGE-DTAB,F,G);  S_PLUS := F*F+G*G-K_SQR;
    REPEAT
      HOUR := HOUR + 2.0*DTAB;
      S_MINUS := S_PLUS;
      FG (    HOUR  ,F,G );   S_0    := F*F+G*G-K_SQR;
      FG ( HOUR+DTAB,F,G );   S_PLUS := F*F+G*G-K_SQR;
      QUAD ( S_MINUS,S_0,S_PLUS, XE,YE, Z[1],Z[2],NZ );
      FOR I:=1 TO NZ DO  TIME[NFOUND+I] := HOUR+DTAB*Z[I];
      NFOUND := NFOUND + NZ;   EVENT := (NFOUND=2);
    UNTIL ( (EVENT) OR (HOUR>=RANGE) ) ;

    (* calculate details of an occultation *)
    IF EVENT THEN
      BEGIN
        MJD_UT_IN  := MJD_CONJ_UT +  TIME[1] / 24.0;
        MJD_UT_OUT := MJD_CONJ_UT +  TIME[2] / 24.0;
        CONTACT ( TIME[1], POS_IN, H_IN, A_IN, B_IN  );
        CONTACT ( TIME[2], POS_OUT,H_OUT,A_OUT,B_OUT );
      END;


  END;

(*-----------------------------------------------------------------------*)
(* DARKNESS: test for civil twilight                                     *)
(*           MODJD:     Modified Julian Date                             *)
(*           LAMBDA:    geographic longitude (>0 west of Greenwich)      *)
(*           CPHI,SPHI: sine and cosine of the geographic latitude       *)
(*-----------------------------------------------------------------------*)
FUNCTION DARKNESS ( MODJD, LAMBDA,CPHI,SPHI: REAL ): BOOLEAN;
  VAR T,RA,DEC,TAU,SIN_HSUN: REAL;
  BEGIN
    T := (MODJD-51544.5)/36525.0;
    MINI_SUN (T,RA,DEC);
    TAU := 15.0 * (LMST(MODJD,LAMBDA) - RA);
    SIN_HSUN  := SPHI*SN(DEC) + CPHI*CS(DEC)*CS(TAU);
    DARKNESS := ( SIN_HSUN < -0.10 );
  END;


(*-----------------------------------------------------------------------*)
(*  EXAMINE:                                                             *)
(*  checks whether an occultation takes place, calculates the            *)
(*  circumstances and prints the results                                 *)
(*                                                                       *)
(*  T1,T2                : search interval in Julian cent. since J2000   *)
(*  RAPOLY,DEPOLY,RPOLY  : Chebyshev approximations of lunar coordinates *)
(*  ETDIFUT              : ET-UT in sec                                  *)
(*  LAMBDA,RCPHI,RSPHI   : geocentric coordinates of the observer        *)
(*  RA_STAR,DE_STAR      : star's coordinates                            *)
(*  NAME                 : star's name                                   *)
(*-----------------------------------------------------------------------*)

PROCEDURE EXAMINE ( T1,T2: REAL; RAPOLY,DEPOLY,RPOLY: TPOLYNOM;
                    ETDIFUT,LAMBDA,RCPHI,RSPHI,RA_STAR,DE_STAR: REAL;
                    NAME: NAME_STRING);

  CONST H_MIN=5.0;  (* minimum altitude above the horizon (deg) *)

  VAR DAY,MONTH,YEAR,H,M,I                               : INTEGER;
      S,HOUR, T_CONJ_ET, MJD_UT_IN, MJD_UT_OUT           : REAL;
      POS_IN,POS_OUT, H_IN,H_OUT, A_IN,A_OUT, B_IN,B_OUT : REAL;
      CONJ, TAKES_PLACE                                  : BOOLEAN;

  BEGIN

    (* test for conjunction in RA and find time of conjunction *)

    CONJUNCT ( T1,T2, RAPOLY,DEPOLY,RPOLY,
               RA_STAR,DE_STAR, CONJ, T_CONJ_ET );

    IF CONJ THEN

      BEGIN

        (* check a possible occultation for the given observing site    *)
        (* and calculate times of contact, altitudes and longitude and  *)
        (* latitude coefficients                                        *)

        SHADOW ( RAPOLY,DEPOLY,RPOLY, T_CONJ_ET, ETDIFUT,
                 LAMBDA,RCPHI,RSPHI, RA_STAR, DE_STAR, TAKES_PLACE,
                 MJD_UT_IN, MJD_UT_OUT, POS_IN,POS_OUT,
                 H_IN,H_OUT, A_IN,A_OUT, B_IN,B_OUT );

        (* print results if the occultation takes place during the *)
        (* night and high enough above the horizon                 *)

        IF TAKES_PLACE THEN
        IF ( (H_IN>H_MIN) OR (H_OUT>H_MIN) ) THEN
        IF DARKNESS ( (MJD_UT_IN+MJD_UT_OUT)/2.0, LAMBDA,RCPHI,RSPHI )
         THEN
          BEGIN

            (* disappearance *)
            CALDAT ( MJD_UT_IN, DAY,MONTH,YEAR,HOUR );   DMS (HOUR,H,M,S);
            WRITE   ( (YEAR MOD 100):3, '/',
                      MONTH:2, '/', DAY:2, H:5, M:3,
                      TRUNC(S+0.5):3, '   D  ', TRUNC(POS_IN+0.5):5,
                      TRUNC(H_IN+0.5):6,  A_IN:8:1,  B_IN:6:1,' ':3    );
            FOR I:=1 TO NAME_LENGTH DO WRITE(NAME[I]); WRITELN;

            (* reapparence *)
            CALDAT ( MJD_UT_OUT, DAY,MONTH,YEAR,HOUR );  DMS (HOUR,H,M,S);
            WRITELN ( (YEAR MOD 100):3, '/',
                      MONTH:2, '/', DAY:2, H:5, M:3,
                      TRUNC(S+0.5):3, '   R  ', TRUNC(POS_OUT+0.5):5,
                      TRUNC(H_OUT+0.5):6, A_OUT:8:1, B_OUT:6:1         );

          END;

      END;

  END;



(*-----------------------------------------------------------------------*)

BEGIN (* main program *)

  (* read search interval and geographic coordinates                     *)

  GET_INPUT ( T_BEGIN,T_END, ETDIFUT, LAMBDA,PHI );

  (* calculate geocentric coordinates of the observer                    *)

  SITE ( PHI, RCPHI,RSPHI );


  (* search occultations in subsequent time intervals                    *)

  T2 := T_BEGIN;

  REPEAT

    T1:=T2;  T2:=T1+T_SEARCH;

    (* approximate lunar coordinates by Chebyshev polynomials            *)

    T_FIT_MOON ( T1-TOVLAP,T2+TOVLAP,MAX_TP_DEG,RAPOLY,DEPOLY,RPOLY );

    (* print header                                                      *)

    HEADER;

    (* open star catalogue file, read epoch and equinox                  *)

    (* RESET ( OCCINP ); *)                           (* Standard Pascal *)
    ASSIGN (OCCINP,'OCCINP.DAT'); RESET(OCCINP);      (* Turbo Pascal    *)

    READLN ( OCCINP, T_EPOCH, T_EQX );
    T_EQX   := ( T_EQX  -2000.0 ) / 100.0;
    T_EPOCH := ( T_EPOCH-2000.0 ) / 100.0;

    (* calculate transformation matrix between the mean equinox of the   *)
    (* star catalog and the true equinox of the search interval center   *)

    TM := (T1+T2)/2.0;
    PN_MATRIX ( T_EQX, TM, PNMAT );

    (* heliocentric velocity of the earth for calculation of aberration  *)

    ABERRAT ( TM, VX,VY,VZ );


    (* loop through list of stars and search for possible occultations   *)

    WHILE NOT EOF(OCCINP) DO
      BEGIN
        (* read new star coordinates      *)
        GETSTAR ( T_EPOCH, TM, RA_STAR,DE_STAR,NAME );
        (* calculate apparent coordinates *)
        APPARENT ( PNMAT, VX,VY,VZ, RA_STAR,DE_STAR );
        (* check for occultation          *)
        EXAMINE ( T1,T2, RAPOLY,DEPOLY,RPOLY,
                  ETDIFUT, LAMBDA,RCPHI,RSPHI, RA_STAR,DE_STAR,NAME );
      END;
    WRITELN;

    CLOSE (OCCINP);  (* Close input file; Turbo Pascal only *)

  UNTIL (T2 >= T_END);

END.

