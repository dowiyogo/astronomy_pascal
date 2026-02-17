(*-----------------------------------------------------------------------*)
(* ELLIP: calculation of position and velocity vector                    *)
(*        for elliptic orbits                                            *)
(*                                                                       *)
(*        M    mean anomaly in degrees       X,Y    position   (in AU)   *)
(*        A    semimajor axis                VX,VY  velocity (in AU/day) *)
(*        ECC  eccentricity                                              *)
(*-----------------------------------------------------------------------*)
PROCEDURE ELLIP(M,A,ECC:REAL;VAR X,Y,VX,VY:REAL);
  CONST KGAUSS = 0.01720209895;
  VAR K,E,C,S,FAC,RHO: REAL;
  BEGIN
    K  := KGAUSS / SQRT(A);
    E  := ECCANOM(M,ECC);   C:=CS(E); S:=SN(E);
    FAC:= SQRT((1.0-ECC)*(1+ECC));  RHO:=1.0-ECC*C;
    X := A*(C-ECC); Y :=A*FAC*S;   VX:=-K*S/RHO;   VY:=K*FAC*C/RHO;
  END;

