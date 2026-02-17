(*-----------------------------------------------------------------------*)
(* HYPERB: calculation of the position and velocity vector               *)
(*         for hyperbolic orbits                                         *)
(*                                                                       *)
(*   T0   time of perihelion passage             X,Y   position          *)
(*   T    time                                   VX,VY velocity          *)
(*   A    semimajor axis (arbitrary sign)                                *)
(*   ECC  eccentricity                                                   *)
(*   (T0,T in julian centuries since J2000)                              *)
(*-----------------------------------------------------------------------*)
PROCEDURE HYPERB(T0,T,A,ECC:REAL;VAR X,Y,VX,VY:REAL);
  CONST KGAUSS = 0.01720209895;
  VAR K,MH,H,EXH,COSHH,SINHH,RHO,FAC: REAL;
  BEGIN
    A   := ABS(A);
    K   := KGAUSS / SQRT(A);
    MH  := K*36525.0*(T-T0)/A;
    H   := HYPANOM(MH,ECC);
    EXH := EXP(H);   COSHH:=0.5*(EXH+1.0/EXH);  SINHH:=0.5*(EXH-1.0/EXH);
    FAC := SQRT((ECC+1.0)*(ECC-1.0));   RHO := ECC*COSHH-1.0;
    X   := A*(ECC-COSHH);  Y  := A*FAC*SINHH;
    VX  :=-K*SINHH/RHO;    VY := K*FAC*COSHH/RHO;
  END;

