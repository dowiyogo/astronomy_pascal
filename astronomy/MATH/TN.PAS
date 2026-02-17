(*-----------------------------------------------------------------------*)
(* TN: tangent function (degrees)                                        *)
(*-----------------------------------------------------------------------*)
FUNCTION TN(X:REAL):REAL;
  CONST RAD=0.0174532925199433;
  VAR XX: REAL;
  BEGIN
    XX:=X*RAD; TN:=SIN(XX)/COS(XX);
  END;
