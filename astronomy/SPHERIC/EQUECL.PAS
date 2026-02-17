(*-----------------------------------------------------------------------*)
(* ECLEQU: Conversion of ecliptic into equatorial coordinates            *)
(*         (T: equinox in Julian centuries since J2000)                  *)
(*-----------------------------------------------------------------------*)
PROCEDURE EQUECL(T:REAL;VAR X,Y,Z:REAL);
  VAR EPS,C,S,V: REAL;
  BEGIN
    EPS:=23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0;
    C:=CS(EPS);  S:=SN(EPS);
    V:=+C*Y+S*Z;  Z:=-S*Y+C*Z;  Y:=V;
  END;
