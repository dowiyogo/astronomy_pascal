(*-----------------------------------------------------------------------*)
(*                                 SUNSET                                *)
(*                 solar and lunar rising and setting times              *)
(*                            version 93/07/01                           *)
(*-----------------------------------------------------------------------*)

PROGRAM SUNSET(INPUT,OUTPUT);

  uses matlib, sunlib, moolib, timlib;

  VAR   ABOVE,RISE,SETT                      : BOOLEAN;
        DAY,MONTH,YEAR, I,IOBJ,NZ            : INTEGER;
        LAMBDA,ZONE,PHI,SPHI,CPHI            : REAL;
        TSTART,DATE,HOUR,HH,UTRISE,UTSET     : REAL;
        Y_MINUS,Y_0,Y_PLUS,ZERO1,ZERO2,XE,YE : REAL;
        SINH0                                : ARRAY[1..3] OF REAL;

(*-----------------------------------------------------------------------*)
(* SIN_ALT: sin(altitude)                                                *)
(*         IOBJ:  1=moon, 2=sun                                          *)
(*-----------------------------------------------------------------------*)
FUNCTION SIN_ALT(IOBJ:INTEGER;MJD0,HOUR,LAMBDA,CPHI,SPHI:REAL): REAL;
  VAR MJD,T,RA,DEC,TAU: REAL;
  BEGIN
    MJD := MJD0 + HOUR/24.0;
    T   := (MJD-51544.5)/36525.0;
    IF (IOBJ=1)
      THEN  MINI_MOON(T,RA,DEC)
      ELSE  MINI_SUN (T,RA,DEC);
    TAU := 15.0 * (LMST(MJD,LAMBDA) - RA);
    SIN_ALT  := SPHI*SN(DEC) + CPHI*CS(DEC)*CS(TAU);
  END;
(*-----------------------------------------------------------------------*)
(* WHM: write time in hours and minutes                                  *)
(*-----------------------------------------------------------------------*)
PROCEDURE WHM(UT:REAL);
  VAR H,M:INTEGER;
  BEGIN
    UT := TRUNC(UT*60.0+0.5)/60.0; (* round to 1 min *)
    H:=TRUNC(UT); M:=TRUNC(60.0*(UT-H)+0.5); WRITE(H:5,':',M:2,' ');
  END;
(*-----------------------------------------------------------------------*)
(* GETINP: read input data                                               *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETINP(VAR DATE,LAMBDA,PHI,ZONE: REAL);
  VAR D,M,Y: INTEGER;
  BEGIN
    WRITELN;
    WRITELN('      SUNSET: solar and lunar rising and setting times ');
    WRITELN('                      version 93/07/01                 ');
    WRITELN('        (c) 1993 Thomas Pfleger, Oliver Montenbruck    ');
    WRITELN;
    WRITELN;
    WRITE  (' First date (yyyy mm dd)               ... '); READLN(Y,M,D);
    WRITELN;
    WRITE  (' Observing site:  longitude (l<0 east) ... '); READLN(LAMBDA);
    WRITE  ('                  latitude             ... '); READLN(PHI);
    WRITE  ('                  local time - UT (h)  ... '); READLN(ZONE);
    WRITELN;
    WRITELN; WRITE('    Date');
    WRITELN('             Moon              Sun             Twilight ');
    WRITELN; WRITE('        ');
    WRITELN('           rise/set          rise/set        end/beginning');
    WRITELN;
    ZONE := ZONE /24.0;
    DATE := MJD(D,M,Y,0)-ZONE;
  END;
(*-----------------------------------------------------------------------*)

BEGIN (* SUNSET *)

  SINH0[1] := SN ( +8.0/60.0); (* moonrise          at h= +8'        *)
  SINH0[2] := SN (-50.0/60.0); (* sunrise           at h=-50'        *)
  SINH0[3] := SN (   -12.0  ); (* nautical twilight at h=-12 degrees *)

  GETINP (TSTART,LAMBDA,PHI,ZONE);

  SPHI := SN(PHI);  CPHI := CS(PHI);

  FOR I:=0 TO 9 DO  (* loop over 10 subsequent days *)

     BEGIN

       DATE := TSTART + I;
       CALDAT(DATE+ZONE,DAY,MONTH,YEAR,HH);
       WRITE(YEAR:5,'/',MONTH:2,'/',DAY:2,'  '); (* print current date *)

       FOR IOBJ := 1 TO 3 DO

         BEGIN

           HOUR := 1.0;

           Y_MINUS := SIN_ALT(IOBJ,DATE,HOUR-1.0,LAMBDA,CPHI,SPHI)
                       - SINH0[IOBJ];

           ABOVE := (Y_MINUS>0.0); RISE := FALSE; SETT := FALSE;

           (* loop over search intervals from [0h-2h] to [22h-24h]  *)
           REPEAT

             Y_0    := SIN_ALT(IOBJ,DATE,HOUR,LAMBDA,CPHI,SPHI)
                         -  SINH0[IOBJ];
             Y_PLUS := SIN_ALT(IOBJ,DATE,HOUR+1.0,LAMBDA,CPHI,SPHI)
                         -  SINH0[IOBJ];

             (* find parabola through three values Y_MINUS,Y_0,Y_PLUS *)
             QUAD ( Y_MINUS,Y_0,Y_PLUS, XE,YE, ZERO1,ZERO2, NZ );

             CASE (NZ) OF
               0: ;
               1: IF (Y_MINUS<0.0)
                    THEN BEGIN UTRISE:=HOUR+ZERO1;  RISE:=TRUE; END
                    ELSE BEGIN UTSET :=HOUR+ZERO1;  SETT:=TRUE; END;
               2: BEGIN
                    IF (YE<0.0)
                      THEN BEGIN UTRISE:=HOUR+ZERO2; UTSET:=HOUR+ZERO1; END
                      ELSE BEGIN UTRISE:=HOUR+ZERO1; UTSET:=HOUR+ZERO2; END;
                    RISE:=TRUE; SETT:=TRUE;
                  END;
              END;

              Y_MINUS := Y_PLUS;     (* prepare for next interval *)
              HOUR := HOUR + 2.0;

           UNTIL ( (HOUR=25.0) OR (RISE AND SETT) );

           IF (RISE OR SETT) (* output *)
             THEN
               BEGIN
                 IF RISE THEN WHM(UTRISE) ELSE WRITE('----- ':9);
                 IF SETT THEN WHM(UTSET)  ELSE WRITE('----- ':9);
               END
             ELSE
               BEGIN
                 IF ABOVE
                   THEN CASE IOBJ OF
                          1,2: WRITE ('   always visible ');
                          3:   WRITE ('    always bright ');
                         END
                   ELSE CASE IOBJ OF
                          1,2: WRITE ('  always invisible');
                          3:   WRITE ('     always dark  ');
                         END;
               END;

         END;

       WRITELN;

     END; (* end of loop over 10 days *)


  WRITELN;  WRITE  (' all times in local time ( = UT ');
  IF ZONE>=0 THEN WRITE('+'); WRITELN(ZONE*24.0:5:1,'h )');

END.


