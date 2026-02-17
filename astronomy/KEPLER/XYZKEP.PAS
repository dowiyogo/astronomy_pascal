(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* XYZKEP: conversion of the state vector into Keplerian elements        *)
(*         for elliptical orbits                                         *)
(*                                                                       *)
(*  X,Y,Z    : Position [AU]                                             *)
(*  VX,VY,VZ : Velocity [AU/d]                                           *)
(*  AX       : Semi-major axis [AU]                                      *)
(*  ECC      : Eccentricity                                              *)
(*  INC      : Inclination [deg]                                         *)
(*  LAN      : Longitude of the ascending node [deg]                     *)
(*  AOP      : Argument of perihelion [deg]                              *)
(*  MA       : Mean anomaly [deg]                                        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE XYZKEP ( X,Y,Z, VX,VY,VZ           :  REAL;
                   VAR AX,ECC,INC,LAN,AOP,MA :  REAL );

  CONST DEG    = 57.29577951308;    (* Conversion from radian to degrees *)
        KGAUSS =  0.01720209895;    (* Gaussian gravitational constant   *)

  VAR   HX,HY,HZ,H, R,V2    : REAL;
        GM, C,S,E2, EA,U,NU : REAL;

  BEGIN

    HX := Y*VZ-Z*VY;                               (* Areal velocity     *)
    HY := Z*VX-X*VZ;
    HZ := X*VY-Y*VX;
    H  := SQRT ( HX*HX + HY*HY + HZ*HZ );

    LAN := ATN2 ( HX, -HY );                       (* Long. ascend. node *)
    INC := ATN2 ( SQRT(HX*HX+HY*HY), HZ );         (* Inclination        *)
    U   := ATN2 ( Z*H, -X*HY+Y*HX );               (* Arg. of latitude   *)

    GM := KGAUSS*KGAUSS;
    R  := SQRT ( X*X + Y*Y + Z*Z );                (* Distance           *)
    V2 := VX*VX + VY*VY + VZ*VZ;                   (* Velocity squared   *)

    AX := 1.0 / (2.0/R-V2/GM);                     (* Semi-major axis    *)

    C  := 1.0-R/AX;                                (* e*cos(E)           *)
    S  := (X*VX+Y*VY+Z*VZ)/(SQRT(AX)*KGAUSS);      (* e*sin(E)           *)

    E2  := C*C + S*S;
    ECC := SQRT ( E2 );                            (* Eccentricity       *)
    EA  := ATN2 ( S, C );                          (* Eccentric anomaly  *)

    MA  := EA - S*DEG;                             (* Mean anomaly       *)

    NU  := ATN2 ( SQRT(1.0-E2)*S, C-E2 );          (* True anomaly       *)
    AOP := U - NU;                                 (* Arg. of perihelion *)

    IF (LAN<0.0) THEN LAN:=LAN+360.0;
    IF (AOP<0.0) THEN AOP:=AOP+360.0;
    IF (MA <0.0) THEN MA :=MA +360.0;

  END;
 
