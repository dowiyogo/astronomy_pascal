(*-----------------------------------------------------------------------*)
(* PN_MATRIX: combined precession and nutation matrix for transformation *)
(*            from mean equinox T0 to true equinox T                     *)
(*            T0,T in Julian cent. since J2000; T=(JD-2451545.0)/36525   *)
(*-----------------------------------------------------------------------*)
PROCEDURE PN_MATRIX ( T0,T:REAL; VAR A: REAL33 );
  BEGIN
    PMATEQU(T0,T,A);                  (* precession matrix T0->T;     *)
    NUTEQU(T,A[1,1],A[2,1],A[3,1]);   (* transform column vectors of  *)
    NUTEQU(T,A[1,2],A[2,2],A[3,2]);   (* matrix A from mean equinox T *)
    NUTEQU(T,A[1,3],A[2,3],A[3,3]);   (* to true equinox T            *)
  END;

