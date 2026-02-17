(*-----------------------------------------------------------------------*)
(* HOREQU: conversion of horizontal to equatorial coordinates            *)
(*   H,AZ : azimuth and altitude (in deg)                                *)
(*   PHI  : geographical latitude (in deg)                               *)
(*   DEC  : declination (-90 deg .. +90 deg)                             *)
(*   TAU  : hour angle (0 deg .. 360 deg)                                *)
(*-----------------------------------------------------------------------*)
PROCEDURE HOREQU (H,AZ,PHI: REAL; VAR DEC,TAU: REAL);
  VAR CS_PHI,SN_PHI, CS_H,SN_H, CS_AZ, X,Y,Z, DUMMY: REAL;
  BEGIN
    CS_PHI:=CS(PHI); SN_PHI:=SN(PHI);
    CS_H  :=CS(H);   SN_H  :=SN(H); CS_AZ:=CS(AZ);
    X := CS_H*SN_PHI*CS_AZ+SN_H*CS_PHI;
    Y := CS_H*SN(AZ);
    Z := SN_H*SN_PHI-CS_H*CS_PHI*CS_AZ;
    POLAR (X,Y,Z, DUMMY,DEC,TAU)
  END; 

