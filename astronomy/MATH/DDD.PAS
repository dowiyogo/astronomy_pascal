(*-----------------------------------------------------------------------*)
(* DDD: conversion of degrees, minutes and seconds into                  *)
(*      degrees and fractions of a degree                                *)
(*-----------------------------------------------------------------------*)
PROCEDURE DDD(D,M:INTEGER;S:REAL;VAR DD:REAL);
 VAR SIGN: REAL;
  BEGIN
    IF ( (D<0) OR (M<0) OR (S<0) ) THEN SIGN:=-1.0 ELSE SIGN:=1.0;
    DD:=SIGN*(ABS(D)+ABS(M)/60.0+ABS(S)/3600.0);
  END;
