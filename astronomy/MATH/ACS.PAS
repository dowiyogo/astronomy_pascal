(*-----------------------------------------------------------------------*)
(* ACS: arccosine function (degrees)                                     *)
(*-----------------------------------------------------------------------*)
FUNCTION ACS(X:REAL):REAL;
  CONST RAD=0.0174532925199433; EPS=1E-7; C=90.0;
  BEGIN
    IF ABS(X)=1.0
     THEN ACS:=C-X*C
     ELSE IF (ABS(X)>EPS) THEN ACS:=C-ARCTAN(X/SQRT((1.0-X)*(1.0+X)))/RAD
                           ELSE ACS:=C-X/RAD;
  END;
