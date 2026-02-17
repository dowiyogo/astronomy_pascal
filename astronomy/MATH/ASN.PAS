(*-----------------------------------------------------------------------*)
(* ASN: arcsine function (degrees)                                       *)
(*-----------------------------------------------------------------------*)
FUNCTION ASN(X:REAL):REAL;
  CONST RAD=0.0174532925199433; EPS=1E-7;
  BEGIN
    IF ABS(X)=1.0
      THEN ASN:=90.0*X
      ELSE IF (ABS(X)>EPS) THEN ASN:=ARCTAN(X/SQRT((1.0-X)*(1.0+X)))/RAD
                           ELSE ASN:=X/RAD;
  END;
