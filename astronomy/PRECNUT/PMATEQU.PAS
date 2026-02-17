(*-----------------------------------------------------------------------*)
(* PMATEQU: calculates the precession matrix A[i,j] for                  *)
(*          transforming equatorial coordinates from equinox T1 to T2    *)
(*          (T=(JD-2451545.0)/36525 )                                    *)
(*-----------------------------------------------------------------------*)
PROCEDURE PMATEQU(T1,T2:REAL;VAR A:REAL33);
  CONST SEC=3600.0;
  VAR DT,ZETA,Z,THETA: REAL;
      C1,S1,C2,S2,C3,S3: REAL;
  BEGIN
   DT:=T2-T1;
    ZETA  :=  ( (2306.2181+(1.39656-0.000139*T1)*T1)+
                ((0.30188-0.000345*T1)+0.017998*DT)*DT )*DT/SEC;
    Z     :=  ZETA + ( (0.79280+0.000411*T1)+0.000205*DT)*DT*DT/SEC;
    THETA :=  ( (2004.3109-(0.85330+0.000217*T1)*T1)-
                ((0.42665+0.000217*T1)+0.041833*DT)*DT )*DT/SEC;
    C1:=CS(Z);  C2:=CS(THETA);  C3:=CS(ZETA);
    S1:=SN(Z);  S2:=SN(THETA);  S3:=SN(ZETA);
    A[1,1]:=-S1*S3+C1*C2*C3; A[1,2]:=-S1*C3-C1*C2*S3; A[1,3]:=-C1*S2;
    A[2,1]:=+C1*S3+S1*C2*C3; A[2,2]:=+C1*C3-S1*C2*S3; A[2,3]:=-S1*S2;
    A[3,1]:=+S2*C3;          A[3,2]:=-S2*S3;          A[3,3]:=+C2;
  END;
