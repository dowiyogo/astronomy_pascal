(*-----------------------------------------------------------------------*)
(* GAUSVEC: calculation of the Gaussian vectors (P,Q,R) from             *)
(*          ecliptic orbital elements:                                   *)
(*          LAN = longitude of the ascending node                        *)
(*          INC = inclination                                            *)
(*          AOP = argument of perihelion                                 *)
(*-----------------------------------------------------------------------*)
PROCEDURE GAUSVEC(LAN,INC,AOP:REAL;VAR PQR:REAL33);
  VAR C1,S1,C2,S2,C3,S3: REAL;
  BEGIN
    C1:=CS(AOP);  C2:=CS(INC);  C3:=CS(LAN);
    S1:=SN(AOP);  S2:=SN(INC);  S3:=SN(LAN);
    PQR[1,1]:=+C1*C3-S1*C2*S3; PQR[1,2]:=-S1*C3-C1*C2*S3; PQR[1,3]:=+S2*S3;
    PQR[2,1]:=+C1*S3+S1*C2*C3; PQR[2,2]:=-S1*S3+C1*C2*C3; PQR[2,3]:=-S2*C3;
    PQR[3,1]:=+S1*S2;          PQR[3,2]:=+C1*S2;          PQR[3,3]:=+C2;
  END;
