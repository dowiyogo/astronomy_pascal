(*-----------------------------------------------------------------------*)
(* EQUHOR: conversion of equatorial into horizontal coordinates          *)
(*   DEC  : declination (-90 deg .. +90 deg)                             *)
(*   TAU  : hour angle (0 deg .. 360 deg)                                *)
(*   PHI  : geographical latitude (in deg)                               *)
(*   H    : altitude (in deg)                                            *)
(*   AZ   : azimuth (0 deg .. 360 deg, counted S->W->N->E->S)            *)
(*-----------------------------------------------------------------------*)
PROCEDURE EQUHOR (DEC,TAU,PHI: REAL; VAR H,AZ: REAL);
  VAR CS_PHI,SN_PHI, CS_DEC,SN_DEC, CS_TAU, X,Y,Z, DUMMY: REAL;
  BEGIN
    CS_PHI:=CS(PHI); SN_PHI:=SN(PHI);
    CS_DEC:=CS(DEC); SN_DEC:=SN(DEC); CS_TAU:=CS(TAU);
    X:=CS_DEC*SN_PHI*CS_TAU - SN_DEC*CS_PHI;
    Y:=CS_DEC*SN(TAU);
    Z:=CS_DEC*CS_PHI*CS_TAU + SN_DEC*SN_PHI;
    POLAR (X,Y,Z, DUMMY,H,AZ)
  END;

