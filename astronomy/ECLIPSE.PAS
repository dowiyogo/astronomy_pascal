(*-----------------------------------------------------------------------*)
(*                                ECLIPSE                                *)
(*             central line and duration of a solar eclipse              *)
(*                          version 93/07/01                             *)
(*-----------------------------------------------------------------------*)
PROGRAM ECLIPSE(INPUT,OUTPUT);

  uses matlib, pnulib, sphlib, sunlib, moolib, timlib;

  CONST MAX_TP_DEG = 8;               (* max. degree of Chebyshev polyn. *)
        H          = 1.14E-6;         (* 1h in jul.cent. (1/(24*36525))  *)

  TYPE  PHASE_TYPE = ( NO_ECLIPSE, PARTIAL,
                       NON_CEN_ANN, NON_CEN_TOT, ANNULAR, TOTAL );

  VAR   T_BEGIN,T_END,T,DT,MJDUT             : REAL;
        ETDIFUT                              : REAL;
        LAMBDA,PHI,T_UMBRA                   : REAL;
        RAM_POLY,DEM_POLY,RM_POLY            : TPOLYNOM;
        RAS_POLY,DES_POLY,RS_POLY            : TPOLYNOM;
        PHASE                                : PHASE_TYPE;


(*-----------------------------------------------------------------------*)
(* GET_INPUT: determines search interval, step size and ET-UT            *)
(* (T_BEGIN,T_END,DT in Julian centuries since J2000 UT; ET-UT in sec)   *)
(*-----------------------------------------------------------------------*)

PROCEDURE GET_INPUT ( VAR T_BEGIN, T_END, DT, ETDIFUT: REAL );

  VAR D,M,Y     : INTEGER;
      UT,T,DTAB : REAL;
      VALID     : BOOLEAN;

  BEGIN

    WRITELN;
    WRITELN ('                     ECLIPSE: solar eclipses            ');
    WRITELN ('                            93/07/01                    ');
    WRITELN ('          (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    WRITELN;
    WRITE   (' Date of New Moon (yyyy mm dd UT): '); READLN(Y,M,D,UT);
    WRITE   (' Output step size (min)          : '); READLN(DTAB);
    DT := (DTAB/1440.0)/36525.0;
    UT := TRUNC ( UT*60.0/DTAB + 0.5 ) * DTAB/60.0 ; (* round to 1 min *)
    T := (MJD(D,M,Y,UT)-51544.5)/36525.0;
    T_BEGIN := T-0.25/36525.0;   T_END := T+0.25/36525.0;
    ETMINUT ( T, ETDIFUT, VALID);
    WRITE   (' Difference ET-UT (sec)          : ');
    IF (VALID) THEN  WRITE ('  (proposal:',TRUNC(ETDIFUT+0.5):4,' sec) ');
    READLN(ETDIFUT);
    WRITELN;
    WRITELN ('    Date       UT         Phi     Lambda   Durat  Phase ');
    WRITELN;
    WRITELN ('               h  m        o  ''      o  ''   min       ');

  END;


(*-----------------------------------------------------------------------*)
(* WRTOUT: formated output                                               *)
(*-----------------------------------------------------------------------*)

PROCEDURE WRTOUT ( MJDUT,LAMBDA,PHI,T_UMBRA: REAL; PHASE: PHASE_TYPE );

  VAR DAY,MONTH,YEAR,H,M: INTEGER;
      HOUR,S            : REAL;

  BEGIN

    CALDAT ( MJDUT,DAY,MONTH,YEAR,HOUR );      (* date  *)
    WRITE  ( YEAR:5,'/',MONTH:2,'/',DAY:2);
    DMS(HOUR+0.5/60.0,H,M,S); WRITE (H:4,M:3); (* time rounded to 1 min *)
    IF ( ORD(PHASE)<ORD(ANNULAR) )
      THEN
        WRITE('       -- --     -- --    ---')
      ELSE
        BEGIN
          DMS(PHI,H,M,S);    WRITE (H:9,M:3);
          DMS(LAMBDA,H,M,S); WRITE (H:7,M:3); WRITE (T_UMBRA:7:1);
        END;
    CASE PHASE OF
      NO_ECLIPSE : WRITE('   ----                  ');
      PARTIAL    : WRITE('   partial               ');
      NON_CEN_ANN: WRITE('   annular (non-central) ');
      NON_CEN_TOT: WRITE('   total (non-central)   ');
      ANNULAR    : WRITE('   annular               ');
      TOTAL      : WRITE('   total                 ');
    END;
    WRITELN;
  END;

(*-----------------------------------------------------------------------*)
(* INTSECT: calculates the intersection of the shadow axis with the      *)
(*          surface of the Earth                                         *)
(*                                                                       *)
(*   RAM,DEM,RM,  equatorial coordinates of Moon and Sun (right asc. RA  *)
(*   RAS,DES,RS:  and declination in deg; distance in Earth radii)       *)
(*   X,Y,Z:       equatorial coord. of the shadow point (in Earth radii) *)
(*   EX,EY,EZ:    unit vector of the shadow axis                         *)
(*   D_UMBRA:     umbra diameter in Earth radii                          *)
(*   PHASE:       phase of the eclipse                                   *)
(*-----------------------------------------------------------------------*)

PROCEDURE INTSECT ( RAM,DEM,RM, RAS,DES,RS: REAL;
                    VAR X,Y,Z, EX,EY,EZ, D_UMBRA: REAL;
                    VAR PHASE: PHASE_TYPE );

  CONST FAC = 0.996633;     (* ratio polar/equatorial Earth radius       *)
        D_M =   0.5450;     (* lunar diameter in units of 1 Earth radius *)
        D_S = 218.25;       (* solar diameter in units of 1 Earth radius *)

  VAR XM,YM,ZM, XS,YS,ZS, XMS,YMS,ZMS, RMS: REAL;
      DELTA, R0, S0, S, D_PENUMBRA        : REAL;

  BEGIN

    CART(RM,DEM,RAM,XM,YM,ZM); ZM:=ZM/FAC; (* solar and lunar coordinat. *)
    CART(RS,DES,RAS,XS,YS,ZS); ZS:=ZS/FAC; (* scale z-coordinate         *)

    XMS:=XM-XS;  YMS:=YM-YS;  ZMS:=ZM-ZS;   (* vector Sun -> Moon,       *)
    RMS := SQRT(XMS*XMS+YMS*YMS+ZMS*ZMS);   (* distance Sun -> Moon      *)
    EX:=XMS/RMS;  EY:=YMS/RMS;  EZ:=ZMS/RMS;(* unit vector Sun -> Moon   *)

    S0 := -( XM*EX + YM*EY + ZM*EZ ); (* dist. Moon -> fundamental plane *)
    DELTA := S0*S0+1.0-XM*XM-YM*YM-ZM*ZM;
    R0 := SQRT(1.0-DELTA);            (* dist.center Earth - shadow axis *)

    D_UMBRA    := (D_S-D_M)*S0/RMS-D_M;     (* diameter of pen-/umbra    *)
    D_PENUMBRA := (D_S+D_M)*S0/RMS+D_M;     (* on the fundamental plane  *)


    (* determine phase and shadow coordinates if required                *)
    IF ( R0 < 1.0 )
      THEN                       (* shadow axis intersects the Earth     *)
        BEGIN                    (* -> total or annular eclipse          *)
          S := S0-SQRT(DELTA);
          D_UMBRA := (D_S-D_M)*(S/RMS)-D_M;   (* umbra diameter on Earth *)
          X:=XM+S*EX;  Y:=YM+S*EY; Z:=ZM+S*EZ;
          Z:=Z*FAC;                           (* rescale z-coordinate    *)
          IF D_UMBRA>0 THEN PHASE:=ANNULAR ELSE PHASE:=TOTAL;
        END
      ELSE
        IF ( R0 < 1.0+0.5*ABS(D_UMBRA) )
          THEN                                (* non-central eclipse     *)
            IF D_UMBRA>0 THEN PHASE:=NON_CEN_ANN ELSE PHASE:=NON_CEN_TOT
          ELSE
            IF ( R0 < 1.0+0.5*D_PENUMBRA)
              THEN PHASE := PARTIAL           (* partial eclipse         *)
              ELSE PHASE := NO_ECLIPSE;       (* no eclipse              *)


  END;


(*-----------------------------------------------------------------------*)
(* CENTRAL: central line, phase and duration of the eclipse              *)
(*                                                                       *)
(*   T_UT:        time in Julian centuries since J2000 UT                *)
(*   ETDIFUT:     difference Ephemeris Time - Universal Time (in sec)    *)
(*   RAM_POLY,DEM_POLY,RM_POLY, RAS_POLY,DES_POLY,RS_POLY:               *)
(*                Chebyshev coefficients for solar and lunar coordinates *)
(*   LAMBDA, PHI: geographic long. and latit. of the shadow center (deg) *)
(*   T_UMBRA:     duration of the total or annular phase (min)           *)
(*   PHASE:       phase of the eclipse                                   *)
(*-----------------------------------------------------------------------*)

PROCEDURE CENTRAL ( T_UT, ETDIFUT            : REAL;
                    RAM_POLY,DEM_POLY,RM_POLY: TPOLYNOM;
                    RAS_POLY,DES_POLY,RS_POLY: TPOLYNOM;
                    VAR LAMBDA,PHI,T_UMBRA   : REAL;
                    VAR PHASE                : PHASE_TYPE );

  CONST AU    = 23454.78;     (* 1AU in earth radii (149597870/6378.14)  *)
        DT    = 0.1;          (* small time interval; dt = 0.1 min       *)
        MPC   = 52596000.0;   (* minutes per Julian century (1440*36525) *)
        OMEGA = 4.3755E-3;    (* angular velocity of the earth (rad/min) *)

  VAR RAM,DEM,RM, RAS,DES,RS, RA,DEC,R, DX,DY,DZ,D, MJDUT  : REAL;
      T,X,Y,Z,EX,EY,EZ,D_UMBRA, XX,YY,ZZ,EXX,EYY,EZZ,DU, W : REAL;
      PH                                                   : PHASE_TYPE;

  (* calculate lunar and solar coordinates from Chebyshev coefficients   *)
  PROCEDURE POSITION ( T: REAL; VAR RAM,DEM,RM, RAS,DES,RS : REAL );
    BEGIN
      RAM:=T_EVAL(RAM_POLY,T); RAS:=T_EVAL(RAS_POLY,T);
      DEM:=T_EVAL(DEM_POLY,T); DES:=T_EVAL(DES_POLY,T);
      RM :=T_EVAL(RM_POLY,T);  RS :=T_EVAL(RS_POLY,T);  RS:=RS*AU;
    END;

  BEGIN

    (* julian centuries since J2000 ET *)
    T := T_UT + ETDIFUT/(86400.0*36525.0);

    (* phase of eclipse and coordinates of the shadow at time T *)
    POSITION ( T, RAM,DEM,RM, RAS,DES,RS );
    INTSECT ( RAM,DEM,RM, RAS,DES,RS, X,Y,Z,EX,EY,EZ, D_UMBRA,PHASE );

    (* for central phase only:  geogr. coord. and duration of totality *)
    IF ( ORD(PHASE) < ORD(ANNULAR) )
      THEN BEGIN LAMBDA:=0.0; PHI:=0.0; T_UMBRA:=0.0; END
      ELSE
        BEGIN
          (* geographic coordinates: *)
          MJDUT := 36525.0*T_UT + 51544.5;
          POLAR ( X,Y,Z, R,DEC,RA );
          PHI    := DEC + 0.1924*SN(2.0*DEC);
          LAMBDA := 15.0*LMST(MJDUT,0.0)-RA;
          IF LAMBDA>+180.0 THEN LAMBDA:=LAMBDA-360.0;
          IF LAMBDA<-180.0 THEN LAMBDA:=LAMBDA+360.0;
          (* duration of totality for this place           *)
          (* (a) shadow coordinates at time T+DT (or T-DT) *)
          POSITION ( T+DT/MPC, RAM,DEM,RM, RAS,DES,RS );  W:=+DT*OMEGA;
          INTSECT ( RAM,DEM,RM, RAS,DES,RS, XX,YY,ZZ,EXX,EYY,EZZ, DU, PH );
          IF (ORD(PH)<ORD(ANNULAR)) THEN
            BEGIN
              POSITION ( T-DT/MPC, RAM,DEM,RM,RAS,DES,RS);  W:=-DT*OMEGA;
              INTSECT (RAM,DEM,RM,RAS,DES,RS, XX,YY,ZZ,EXX,EYY,EZZ, DU,PH);
            END;
          (* (b) displacement DX,DY,DZ of the shadow on Earth    *)
          (*     and fraction D perpendicular to the shadow axis *)
          DX := XX-X+W*Y;   DY := YY-Y-W*X;   DZ := ZZ-Z;
          D  := SQRT( DX*DX+DY*DY+DZ*DZ -
                      (DX*EX+DY*EY+DZ*EZ)*(DX*EX+DY*EY+DZ*EZ) );
          T_UMBRA := DT * ABS(D_UMBRA) / D;
        END;

  END;

(*-----------------------------------------------------------------------*)

BEGIN (* main program *)


  (* read search interval *)

  GET_INPUT ( T_BEGIN,T_END, DT, ETDIFUT );


  (* Chebyshev approximations *)

  T_FIT_MOON (T_BEGIN-H,T_END+H,8,RAM_POLY,DEM_POLY,RM_POLY);
  T_FIT_SUN  (T_BEGIN-H,T_END+H,3,RAS_POLY,DES_POLY,RS_POLY);


  (* calculate phase and central line of the eclipse *)

  T := T_BEGIN;

  REPEAT

    CENTRAL ( T,ETDIFUT,
              RAM_POLY,DEM_POLY,RM_POLY, RAS_POLY,DES_POLY,RS_POLY,
              LAMBDA, PHI, T_UMBRA, PHASE  );

    IF PHASE<>NO_ECLIPSE THEN
      BEGIN
        MJDUT := 36525.0*T + 51544.5;
        WRTOUT ( MJDUT, LAMBDA, PHI, T_UMBRA, PHASE );
      END;

    T := T+DT;

  UNTIL (T > T_END);

END.
