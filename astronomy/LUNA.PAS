(*-----------------------------------------------------------------------*)
(*                                 LUNA                                  *)
(*                           lunar ephemeris                             *)
(*                           version 93/07/01                            *)
(*-----------------------------------------------------------------------*)

PROGRAM LUNA(INPUT,OUTPUT);

  uses matlib, pnulib, sphlib, moolib, timlib;

  CONST T_OVERLAP  = 3.42E-6;        (* 3h  in Julian centuries *)
        T_DEVELOP  = 2.737850787E-4; (* 10d in Julian centuries *)

  VAR RA,DE,R,PAR, MODJD,HOUR : REAL;
      T,DT,T_START,T_END,TA,TB: REAL;
      DAY,MONTH,YEAR,NLINE    : INTEGER;
      RA_POLY,DE_POLY,R_POLY  : TPOLYNOM;


(*-----------------------------------------------------------------------*)
(* GETEPH: reads the period of time for the ephemeris                    *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETEPH(VAR T1,DT,T2:REAL);
  VAR YEAR,MONTH,DAY: INTEGER;
      HOUR          : REAL;
  BEGIN
    WRITELN;
    WRITELN('                          LUNA: lunar ephemeris            ');
    WRITELN('                            Version 93/07/01               ');
    WRITELN('             (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    WRITELN;
    WRITELN(' Begin and end of the ephemeris: ');
    WRITELN;
    WRITE  ('  first date (yyyy mm dd hh.hhh)  ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  last date  (yyyy mm dd hh.hhh)  ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  step size (dd hh.hh)            ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
  END;
(*-----------------------------------------------------------------------*)
(* WRTLBRP: formatted output                                             *)
(*-----------------------------------------------------------------------*)
PROCEDURE WRTLBRP (L,B,R,P:REAL);
  VAR H,M: INTEGER;
      S  : REAL;
  BEGIN
    DMS(L,H,M,S);  WRITE (H:5,M:3,S:5:1);
    DMS(B,H,M,S);  WRITE (H:5,M:3,S:5:1);  WRITE (R:10:3);
    DMS(P,H,M,S);  IF (H>0) THEN M:=M+60;  WRITELN (M:6,S:6:2);
  END;
(*-----------------------------------------------------------------------*)

BEGIN   (* main program *)

  GETEPH(T_START,DT,T_END);  (* read desired dates *)

  WRITELN;
  WRITE   ('    Date      ET       RA           Dec      Distance  ');
  WRITELN (' Parallax ');
  WRITE   ('              h      h  m  s      o  ''  "  Earth radii');
  WRITELN ('    ''  "   ');

  T := T_START;   TB := T_START;   NLINE := 0;

  WHILE (T<=T_END) DO

    BEGIN

      IF (T>TB-T_OVERLAP) THEN       (* new expansion of the coordinates *)
        BEGIN
          TA := T-T_OVERLAP;  TB := T+T_DEVELOP+T_OVERLAP;
          T_FIT_MOON (TA,TB,MAX_TP_DEG,RA_POLY,DE_POLY,R_POLY);
        END;

      (* date *)
      MODJD := T*36525.0 + 51544.5;  CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);
      WRITE  (YEAR:5,'/',MONTH:2,'/',DAY:2,HOUR:5:1);

      (* coordinates *)
      RA  := T_EVAL(RA_POLY,T)/15.0;  IF RA<0.0 THEN RA := RA + 24.0;
      DE  := T_EVAL(DE_POLY,T);
      R   := T_EVAL(R_POLY, T);
      PAR := ASN(1.0/R);

      (* print coordinates *)
      WRTLBRP (RA,DE,R,PAR);
      NLINE := NLINE + 1;
      IF (NLINE MOD 5) = 0 THEN WRITELN;

      T := T + DT;

    END;

END.

