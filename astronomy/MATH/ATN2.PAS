(*-----------------------------------------------------------------------*)
(* ATN2: arctangent of y/x for two arguments                             *)
(*       (correct quadrant; -180 deg <= ATN2 <= +180 deg)                *)
(*-----------------------------------------------------------------------*)
FUNCTION ATN2(Y,X:REAL):REAL;
  CONST RAD=0.0174532925199433;
  VAR   AX,AY,PHI: REAL;
  BEGIN
    IF (X=0.0) AND (Y=0.0) 
      THEN ATN2:=0.0
      ELSE 
        BEGIN
          AX:=ABS(X); AY:=ABS(Y); 
          IF (AX>AY) 
            THEN PHI:=ARCTAN(AY/AX)/RAD 
            ELSE PHI:=90.0-ARCTAN(AX/AY)/RAD;
          IF (X<0.0) THEN PHI:=180.0-PHI;
          IF (Y<0.0) THEN PHI:=-PHI;
          ATN2:=PHI;
        END;
  END;
