(*-----------------------------------------------------------------------*)
(* DMS: conversion of degrees and fractions of a degree                  *)
(*      into degrees, minutes and seconds                                *)
(*-----------------------------------------------------------------------*)
PROCEDURE DMS(DD:REAL;VAR D,M:INTEGER;VAR S:REAL);
  VAR D1:REAL; 
  BEGIN
    D1:=ABS(DD);  D:=TRUNC(D1);  
    D1:=(D1-D)*60.0;  M:=TRUNC(D1);  S:=(D1-M)*60.0;
    IF (DD<0) THEN 
      IF (D<>0) THEN D:=-D ELSE IF (M<>0) THEN M:=-M ELSE S:=-S;
  END;
