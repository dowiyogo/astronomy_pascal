(*----------------------------------------------------------------------*)
(* CALDAT: Finds the civil calendar date for a given value              *)
(*         of the Modified Julian Date (MJD).                           *)
(*         Julian calendar is used up to 1582 October 4,                *)
(*         Gregorian calendar is used from 1582 October 15 onwards.     *)
(*----------------------------------------------------------------------*)
PROCEDURE CALDAT(MJD:REAL; VAR DAY,MONTH,YEAR:INTEGER;VAR HOUR:REAL);
  VAR B,D,F     : INTEGER;
      JD,JD0,C,E: REAL;
  BEGIN
    JD  := MJD + 2400000.5;
    (* JD0 := TRUNC(JD+0.5);      *)    (* Standard Pascal  *)
       JD0 := INT(JD+0.5);              (* TURBO Pascal     *)
    (* JD0 := LONG_TRUNC(JD+0.5); *)    (* ST Pascal plus   *)
    IF (JD0<2299161.0)                            (* calendar:    *)
      THEN BEGIN C:=JD0+1524.0 END                (* -> Julian    *)
      ELSE BEGIN                                  (* -> Gregorian *)
             B:=TRUNC((JD0-1867216.25)/36524.25);
             C:=JD0+(B-TRUNC(B/4))+1525.0
           END;
    D    := TRUNC((C-122.1)/365.25);          E     := 365.0*D+TRUNC(D/4);
    F    := TRUNC((C-E)/30.6001);
    DAY  := TRUNC(C-E+0.5)-TRUNC(30.6001*F);  MONTH := F-1-12*TRUNC(F/14);
    YEAR := D-4715-TRUNC((7+MONTH)/10);       HOUR  := 24.0*(JD+0.5-JD0);
  END;

