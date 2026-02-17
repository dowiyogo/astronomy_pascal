(*-----------------------------------------------------------------------*)
(* CROSS: cross product of two vectors                                   *)
(*-----------------------------------------------------------------------*)
PROCEDURE CROSS(A,B:VECTOR;VAR C:VECTOR);
  BEGIN
    C[X] := A[Y]*B[Z]-A[Z]*B[Y];
    C[Y] := A[Z]*B[X]-A[X]*B[Z];
    C[Z] := A[X]*B[Y]-A[Y]*B[X];
  END;

