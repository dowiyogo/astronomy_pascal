(*-----------------------------------------------------------------------*)
(* ORBECL: transformation of coordinates XX,YY referred to the           *)
(*         orbital plane into ecliptic coordinates X,Y,Z using           *)
(*         Gaussian vectors PQR                                          *)
(*-----------------------------------------------------------------------*)
PROCEDURE ORBECL(XX,YY:REAL;PQR:REAL33;VAR X,Y,Z:REAL);
  BEGIN
    X:=PQR[1,1]*XX+PQR[1,2]*YY;
    Y:=PQR[2,1]*XX+PQR[2,2]*YY;
    Z:=PQR[3,1]*XX+PQR[3,2]*YY;
  END;
