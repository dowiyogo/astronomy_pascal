(*-----------------------------------------------------------------------*)
(* PMATECL: calculates the precession matrix A[i,j] for                  *)
(*          transforming ecliptic coordinates from equinox T1 to T2      *)
(*          ( T=(JD-2451545.0)/36525 )                                   *)
(*-----------------------------------------------------------------------*)
PROCEDURE PMATECL(T1,T2:REAL;VAR A: REAL33);
  CONST SEC=3600.0;
  VAR DT,PPI,PI,PA: REAL;
      C1,S1,C2,S2,C3,S3: REAL;
  BEGIN
    DT:=T2-T1;
    PPI := 174.876383889 +( ((3289.4789+0.60622*T1)*T1) +
              ((-869.8089-0.50491*T1) + 0.03536*DT)*DT )/SEC;
    PI  := ( (47.0029-(0.06603-0.000598*T1)*T1)+
             ((-0.03302+0.000598*T1)+0.000060*DT)*DT )*DT/SEC;
    PA  := ( (5029.0966+(2.22226-0.000042*T1)*T1)+
             ((1.11113-0.000042*T1)-0.000006*DT)*DT )*DT/SEC;
    C1:=CS(PPI+PA);  C2:=CS(PI);  C3:=CS(PPI);
    S1:=SN(PPI+PA);  S2:=SN(PI);  S3:=SN(PPI);
    A[1,1]:=+C1*C3+S1*C2*S3; A[1,2]:=+C1*S3-S1*C2*C3; A[1,3]:=-S1*S2;
    A[2,1]:=+S1*C3-C1*C2*S3; A[2,2]:=+S1*S3+C1*C2*C3; A[2,3]:=+C1*S2;
    A[3,1]:=+S2*S3;          A[3,2]:=-S2*C3;          A[3,3]:=+C2;
  END;
